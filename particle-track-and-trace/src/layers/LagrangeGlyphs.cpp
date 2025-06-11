#include "LagrangeGlyphs.h"
#include "../commands/SpawnPointCallback.h"
#include <vtkActor2D.h>
#include <vtkGlyph2D.h>
#include <vtkLookupTable.h>
#include <vtkNamedColors.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper2D.h>
#include <vtkProperty.h>
#include <vtkProperty2D.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkInteractorStyle.h>
#include <vtkInteractorStyleUser.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkCamera.h>
#include <fstream>
#include <sstream>
#include <iostream>

#include "../CartographicTransformation.h"
#include "../advection/interpolate.h"
#include "../advection/kernel/SnapBoundaryConditionKernel.h"
#include "../advection/kernel/FreeSlipBoundaryConditionKernel.h"
#include "../advection/kernel/PartialSlipBoundaryConditionKernel.h"
#include "../advection/kernel/ParcelsBoundaryConditionKernel.h"
#include <cmath>
#include <algorithm>
#include <limits>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

vtkSmartPointer<SpawnPointCallback> LagrangeGlyphs::createSpawnPointCallback()
{
    auto newPointCallBack = vtkSmartPointer<SpawnPointCallback>::New();
    newPointCallBack->setData(this->data);
    newPointCallBack->setPoints(this->points);
    newPointCallBack->setRen(this->renderer);
    newPointCallBack->setUVGrid(this->uvGrid);
    newPointCallBack->setBeached(this->particlesBeached);
    return newPointCallBack;
}

LagrangeGlyphs::LagrangeGlyphs(std::shared_ptr<UVGrid> grid, std::unique_ptr<AdvectionKernel> advectionKernel, const std::string &spawnFile) : uvGrid{std::move(grid)}, advector{std::move(advectionKernel)}
{
    this->data = vtkSmartPointer<vtkPolyData>::New();
    this->data->SetPoints(this->points);

    this->particlesBeached = vtkSmartPointer<vtkIntArray>::New();
    this->particlesBeached->SetName("particlesBeached");
    this->particlesBeached->SetNumberOfComponents(1);
    this->data->GetPointData()->SetScalars(this->particlesBeached);

    vtkSmartPointer<vtkTransformFilter> transformFilter = createCartographicTransformFilter(*uvGrid);
    transformFilter->SetInputData(data);

    circleSource->SetGlyphTypeToCircle();
    circleSource->SetScale(0.01);
    circleSource->Update();

    vtkNew<vtkGlyph2D> glyph2D;
    glyph2D->SetSourceConnection(circleSource->GetOutputPort());
    glyph2D->SetInputConnection(transformFilter->GetOutputPort());
    glyph2D->SetScaleModeToDataScalingOff();
    glyph2D->Update();

    // Create lookup table for coloring particles
    vtkNew<vtkLookupTable> lut;
    lut->SetNumberOfTableValues(2);       // Two states: not beached and beached
    lut->SetTableValue(0, 0.0, 0.0, 1.0); // Blue for not beached
    lut->SetTableValue(1, 1.0, 0.0, 0.0); // Red for beached
    lut->Build();

    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputConnection(glyph2D->GetOutputPort());
    mapper->SetLookupTable(lut);
    mapper->SetScalarRange(0, 1);
    mapper->Update();

    actor->SetMapper(mapper);
    this->renderer->AddActor(actor);

    // Load and spawn particles from CSV file
    spawnParticlesFromFile(spawnFile);

    // Initialize dualConditionResidenceTimes
    dualConditionResidenceTimes.resize(this->points->GetNumberOfPoints(), 0);
}

void LagrangeGlyphs::spawnParticlesFromFile(const std::string &filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open spawn locations file: " << filename << std::endl;
        return;
    }

    std::string line;
    // Skip header line if it exists
    std::getline(file, line);

    while (std::getline(file, line))
    {
        std::stringstream ss(line);
        std::string lat_str, lon_str;

        if (std::getline(ss, lat_str, ',') && std::getline(ss, lon_str, ','))
        {
            try
            {
                double lat = std::stod(lat_str);
                double lon = std::stod(lon_str);

                // Only spawn if within grid boundaries
                if (lon >= uvGrid->lonMin() && lon <= uvGrid->lonMax() &&
                    lat >= uvGrid->latMin() && lat <= uvGrid->latMax())
                {
                    this->points->InsertNextPoint(lon, lat, 0);
                    this->particlesBeached->InsertNextValue(0);
                    this->coastalResidenceTimes.push_back(0);
                }
            }
            catch (const std::exception &e)
            {
                std::cerr << "Error parsing line: " << line << std::endl;
                continue;
            }
        }
    }

    if (this->points->GetNumberOfPoints() > 0)
    {
        this->points->Modified();
        this->particlesBeached->Modified();
        std::cout << "Spawned " << this->points->GetNumberOfPoints() << " particles from file." << std::endl;
    }
}

void LagrangeGlyphs::updateData(int t)
{
    const int SUPERSAMPLINGRATE = 4;
    double point[3], oldX, oldY;
    bool modifiedData = false;
    int dtTotal = t - this->lastT;
    if (dtTotal < 0)
        dtTotal = t;

    // iterate over every point.
    for (vtkIdType n = 0; n < this->points->GetNumberOfPoints(); n++)
    {
        int beachedFor = this->particlesBeached->GetValue(n);
        // first check: only update non-beached particles.
        if (beachedFor < this->beachedAtNumberOfTimes)
        {
            this->points->GetPoint(n, point);
            // second check: only update points within our grid's boundary.
            if (point[0] <= this->uvGrid->lonMin() || point[0] >= this->uvGrid->lonMax() ||
                point[1] <= this->uvGrid->latMin() || point[1] >= this->uvGrid->latMax())
            {
                // sets any particle out of bounds to be beached - so it gets assigned the right colour in the lookup table.
                this->particlesBeached->SetValue(n, this->beachedAtNumberOfTimes);
                continue;
            }

            oldX = point[0];
            oldY = point[1];

            // supersampling
            for (int i = 0; i < SUPERSAMPLINGRATE; i++)
            {
                int dt = (t - this->lastT) / SUPERSAMPLINGRATE;
                if (dt < 0)
                {
                    // TODO: This is a hack for when the t wraps around,
                    // there is probably a more elegant way of dealing with this whole thing
                    // that involves having two separate DTs.
                    // One for the "render" step time, and one for the computation step time.
                    dt = t;
                }
                std::tie(point[1], point[0]) = advector->advect(t, point[1], point[0], dt);

                // Track particle if enabled
                if (trackingEnabled && n == trackedParticleIndex && i == SUPERSAMPLINGRATE - 1)
                {
                    trackedPositions.push_back({point[1], point[0]});

                    // Calculate velocity
                    auto vel = bilinearinterpolate(*uvGrid, t, point[1], point[0]);
                    trackedVelocities.push_back({vel.u, vel.v});

                    // Use GEBCO-based distance to shore
                    double shoreDist = uvGrid->getShoreDistance(point[1], point[0]);
                    trackedDistancesToShore.push_back(shoreDist);
                }

                // Track all particles if enabled
                if (trackingAllEnabled && i == SUPERSAMPLINGRATE - 1)
                {
                    allParticlePositions[n].push_back({point[1], point[0]});

                    // Calculate velocity
                    auto vel = bilinearinterpolate(*uvGrid, t, point[1], point[0]);
                    allParticleVelocities[n].push_back({vel.u, vel.v});

                    // Use GEBCO-based distance to shore
                    double shoreDist = uvGrid->getShoreDistance(point[1], point[0]);
                    allParticleDistancesToShore[n].push_back(shoreDist);
                }
            }

            // Apply beaching logic based on selected type
            switch (beachingType) {
                case BeachingType::VelocityBased:
                    if (isNearestNeighbourZero(*uvGrid, t, point[1], point[0]))
                    {
                        this->particlesBeached->SetValue(n, beachedFor + 1);
                    }
                    else
                    {
                        this->particlesBeached->SetValue(n, std::max(beachedFor - 1, 0));
                        this->points->SetPoint(n, point);
                        modifiedData = true;
                    }
                    break;

                case BeachingType::DistanceBased:
                    {
                        double shoreDist = uvGrid->getShoreDistance(point[1], point[0]);
                        if (shoreDist < 1e-4)
                        {
                            this->particlesBeached->SetValue(n, this->beachedAtNumberOfTimes);
                        }
                        else
                        {
                            this->particlesBeached->SetValue(n, 0);
                            this->points->SetPoint(n, point);
                            modifiedData = true;
                        }
                    }
                    break;

                case BeachingType::DirectionalBased:
                    {
                        const double bufferDistance = 2000.0;     // 2 km bufferzone
                        double shoreDist = uvGrid->getShoreDistance(point[1], point[0]);
                        bool inBuffer = shoreDist < bufferDistance;

                        if (inBuffer && enableDirectionalCheck)
                        {
                            // Calculate shore gradient
                            const double delta = AdvectionKernel::metreToDegrees(100);
                            double gradLat = (uvGrid->getShoreDistance(oldY + delta, oldX) -
                                            uvGrid->getShoreDistance(oldY - delta, oldX)) / (2 * delta);
                            double gradLon = (uvGrid->getShoreDistance(oldY, oldX + delta) -
                                            uvGrid->getShoreDistance(oldY, oldX - delta)) / (2 * delta);

                            // Determine movement direction
                            auto vel = bilinearinterpolate(*uvGrid, t, oldY, oldX);
                            double dot = vel.u * gradLon + vel.v * gradLat;

                            if (dot < 0  && coastalResidenceTimes[n] >= coastalTimeThreshold) // Moving towards shore
                            {
                                this->particlesBeached->SetValue(n, this->beachedAtNumberOfTimes);
                            }
                            else
                            {
                                if (dot >= 0)
                                    coastalResidenceTimes[n] = 0;
                                this->particlesBeached->SetValue(n, 0);
                                this->points->SetPoint(n, point);
                                modifiedData = true;
                            }
                        }
                        else
                        {
                            if (!inBuffer)
                                coastalResidenceTimes[n] = 0;
                            this->particlesBeached->SetValue(n, 0);
                            this->points->SetPoint(n, point);
                            modifiedData = true;
                        }
                    }
                    break;

                case BeachingType::None:
                    // Only beach if outside grid
                    if (point[0] <= uvGrid->lonMin() || point[0] >= uvGrid->lonMax() ||
                        point[1] <= uvGrid->latMin() || point[1] >= uvGrid->latMax())
                    {
                        this->particlesBeached->SetValue(n, this->beachedAtNumberOfTimes);
                    }
                    else
                    {
                        this->particlesBeached->SetValue(n, 0);
                        this->points->SetPoint(n, point);
                        modifiedData = true;
                    }
                    break;
            }

            bool inCoastal = uvGrid->isNearShore(point[1], point[0], 5000.0);
            if (inCoastal)
            {
                if (coastalResidenceTimes.size() > static_cast<size_t>(n))
                    coastalResidenceTimes[n] += dtTotal;
            }
            else if (coastalResidenceTimes.size() > static_cast<size_t>(n))
            {
                coastalResidenceTimes[n] = 0;
            }
        }
    }
    if (modifiedData)
        this->points->Modified();
    this->lastT = t;
}

void LagrangeGlyphs::addObservers(vtkSmartPointer<vtkRenderWindowInteractor> interactor)
{
    auto newPointCallBack = createSpawnPointCallback();
    interactor->AddObserver(vtkCommand::LeftButtonPressEvent, newPointCallBack);
    interactor->AddObserver(vtkCommand::LeftButtonReleaseEvent, newPointCallBack);
    interactor->AddObserver(vtkCommand::MouseMoveEvent, newPointCallBack);
    interactor->AddObserver(vtkCommand::RightButtonPressEvent, newPointCallBack);
}

vtkSmartPointer<vtkPoints> LagrangeGlyphs::getPoints()
{
    return points;
}

vtkSmartPointer<vtkIntArray> LagrangeGlyphs::getBeached()
{
    return particlesBeached;
}

void LagrangeGlyphs::setColour(int red, int green, int blue)
{
    actor->GetProperty()->SetColor(red / 255.0, green / 255.0, blue / 255.0);
}

void LagrangeGlyphs::setToDiamond()
{
    circleSource->SetGlyphTypeToDiamond();
}

void LagrangeGlyphs::handleGameOver()
{
    points->Reset();
    particlesBeached->Reset();
}

void LagrangeGlyphs::startTracking(size_t particleIndex)
{
    trackingEnabled = true;
    trackedParticleIndex = particleIndex;
    trackedPositions.clear();
    trackedVelocities.clear();
    trackedDistancesToShore.clear();
}

void LagrangeGlyphs::stopTracking()
{
    trackingEnabled = false;
}

void LagrangeGlyphs::setEnableDirectionalCheck(bool enabled) {
    enableDirectionalCheck = enabled;
}

void LagrangeGlyphs::setCoastalTimeThreshold(double hours) {
    coastalTimeThreshold = static_cast<int>(hours * 3600.0);
}

void LagrangeGlyphs::printTrackedParticleInfo(const std::string &outputFilename) const
{
    if (!trackingEnabled || trackedPositions.empty())
    {
        std::cout << "No particle is being tracked or no data available." << std::endl;
        return;
    }

    // Save to file
    std::string outputPath = "C:/Users/wesle/Documents/Universiteit/2024_2025/Thesis/Opdracht/interactive-track-and-trace/particle-track-and-trace/" + outputFilename;
    std::ofstream outFile(outputPath);
    if (!outFile.is_open())
    {
        std::cerr << "Error: Could not open file for writing: " << outputPath << std::endl;
        return;
    }

    try
    {
        // Write header
        outFile << "Step,Latitude,Longitude,VelocityU,VelocityV,DistanceToShore,Beached\n";

        // Write data
        for (size_t i = 0; i < trackedPositions.size(); ++i)
        {
            const auto &pos = trackedPositions[i];
            const auto &vel = trackedVelocities[i];
            double shoreDist = trackedDistancesToShore[i];

            // Calculate if particle is beached using the exact logic from DualBoundaryConditionKernel
            bool isBeached = false;
            const double buffer = 5000.0;     // metres
            const double beachThresh = 200.0; // metres

            if (uvGrid->isNearShore(pos.first, pos.second, buffer))
            {
                // Calculate shore gradient
                const double delta = AdvectionKernel::metreToDegrees(100);
                double gradLat = (uvGrid->getShoreDistance(pos.first + delta, pos.second) -
                                  uvGrid->getShoreDistance(pos.first - delta, pos.second)) /
                                 (2 * delta);
                double gradLon = (uvGrid->getShoreDistance(pos.first, pos.second + delta) -
                                  uvGrid->getShoreDistance(pos.first, pos.second - delta)) /
                                 (2 * delta);

                // Calculate dot product
                double dot = vel.first * gradLon + vel.second * gradLat;

                // A particle is beached if:
                // 1. It is moving towards shore (dot < 0)
                // 2. It is within 200m of shore
                isBeached = (dot < 0) && (shoreDist < beachThresh);
            }

            outFile << i << ","
                    << pos.first << ","
                    << pos.second << ","
                    << vel.first << ","
                    << vel.second << ","
                    << shoreDist << ","
                    << (isBeached ? "Yes" : "No") << "\n";
        }

        outFile.close();
        std::cout << "\nTrajectory data saved to: " << outputPath << std::endl;
        std::cout << "Run 'python plot_trajectory.py' to visualize the particle path on a map." << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error writing to file: " << e.what() << std::endl;
        if (outFile.is_open())
        {
            outFile.close();
        }
    }
}

void LagrangeGlyphs::startTrackingAll()
{
    trackingAllEnabled = true;
    allParticlePositions.clear();
    allParticleVelocities.clear();
    allParticleDistancesToShore.clear();

    // Initialize vectors for each particle
    size_t numParticles = points->GetNumberOfPoints();
    allParticlePositions.resize(numParticles);
    allParticleVelocities.resize(numParticles);
    allParticleDistancesToShore.resize(numParticles);
}

void LagrangeGlyphs::stopTrackingAll()
{
    trackingAllEnabled = false;
}

void LagrangeGlyphs::printAllParticlesInfo(const std::string &outputFilename) const
{
    if (!trackingAllEnabled || allParticlePositions.empty())
    {
        std::cout << "No particles are being tracked or no data available." << std::endl;
        return;
    }

    // Save to file
    std::string outputPath = "C:/Users/wesle/Documents/Universiteit/2024_2025/Thesis/Opdracht/interactive-track-and-trace/particle-track-and-trace/" + outputFilename;
    std::ofstream outFile(outputPath);
    if (!outFile.is_open())
    {
        std::cerr << "Error: Could not open file for writing: " << outputPath << std::endl;
        return;
    }

    try
    {
        // Write header
        outFile << "ParticleID,Step,Latitude,Longitude,VelocityU,VelocityV,DistanceToShore\n";

        // Write data for each particle
        for (size_t particleId = 0; particleId < allParticlePositions.size(); ++particleId)
        {
            const auto &positions = allParticlePositions[particleId];
            const auto &velocities = allParticleVelocities[particleId];
            const auto &distances = allParticleDistancesToShore[particleId];

            for (size_t step = 0; step < positions.size(); ++step)
            {
                const auto &pos = positions[step];
                const auto &vel = velocities[step];
                outFile << particleId << ","
                        << step << ","
                        << pos.first << ","
                        << pos.second << ","
                        << vel.first << ","
                        << vel.second << ","
                        << distances[step] << "\n";
            }
        }

        outFile.close();
        std::cout << "\nAll particles trajectory data saved to: " << outputPath << std::endl;
        std::cout << "Run 'python plot_trajectory.py' to visualize the particle paths on a map." << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error writing to file: " << e.what() << std::endl;
        if (outFile.is_open())
        {
            outFile.close();
        }
    }
}
