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
#include <cmath>
#include <algorithm>
#include <limits>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

vtkSmartPointer<SpawnPointCallback> LagrangeGlyphs::createSpawnPointCallback() {
  auto newPointCallBack = vtkSmartPointer<SpawnPointCallback>::New();
  newPointCallBack->setData(this->data);
  newPointCallBack->setPoints(this->points);
  newPointCallBack->setRen(this->renderer);
  newPointCallBack->setUVGrid(this->uvGrid);
  newPointCallBack->setBeached(this->particlesBeached);
  return newPointCallBack;
}

LagrangeGlyphs::LagrangeGlyphs(std::shared_ptr<UVGrid> grid, std::unique_ptr<AdvectionKernel> advectionKernel, const std::string& spawnFile) :
        uvGrid{std::move(grid)}, advector{std::move(advectionKernel)} {

  //Which advection are we using.
  if (dynamic_cast<SnapBoundaryConditionKernel*>(advector.get())) {
      this->boundaryType = BoundaryType::Snap;
  } else if (dynamic_cast<FreeSlipBoundaryConditionKernel*>(advector.get())) {
      this->boundaryType = BoundaryType::FreeSlip;
  } else if (dynamic_cast<PartialSlipBoundaryConditionKernel*>(advector.get())) {
      this->boundaryType = BoundaryType::PartialSlip;
  } else {
      throw std::runtime_error("Unknown boundary kernel type");
  }

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
  lut->SetNumberOfTableValues(2);  // Two states: not beached and beached
  lut->SetTableValue(0, 0.0, 0.0, 1.0);  // Blue for not beached
  lut->SetTableValue(1, 1.0, 0.0, 0.0);  // Red for beached
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
}

void LagrangeGlyphs::spawnParticlesFromFile(const std::string& filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Error: Could not open spawn locations file: " << filename << std::endl;
    return;
  }

  std::string line;
  // Skip header line if it exists
  std::getline(file, line);

  while (std::getline(file, line)) {
    std::stringstream ss(line);
    std::string lat_str, lon_str;
    
    if (std::getline(ss, lat_str, ',') && std::getline(ss, lon_str, ',')) {
      try {
        double lat = std::stod(lat_str);
        double lon = std::stod(lon_str);
        
        // Only spawn if within grid boundaries
        if (lon >= uvGrid->lonMin() && lon <= uvGrid->lonMax() &&
            lat >= uvGrid->latMin() && lat <= uvGrid->latMax()) {
          this->points->InsertNextPoint(lon, lat, 0);
          this->particlesBeached->InsertNextValue(0);
          this->initialSpawnPositions.push_back({lat, lon}); // Store the initial position
        }
      } catch (const std::exception& e) {
        std::cerr << "Error parsing line: " << line << std::endl;
        continue;
      }
    }
  }

  if (this->points->GetNumberOfPoints() > 0) {
    this->points->Modified();
    this->particlesBeached->Modified();
    std::cout << "Spawned " << this->points->GetNumberOfPoints() << " particles from file." << std::endl;
    
    // The initialization of allParticlePositions is now handled in startTrackingAll
  }
}

void LagrangeGlyphs::updateData(int t) {
  const int SUPERSAMPLINGRATE = 4;
  double point[3], oldX, oldY;
  bool modifiedData = false;

  // iterate over every point.
  for (vtkIdType n = 0; n < this->points->GetNumberOfPoints(); n++) {
    int beachedFor = this->particlesBeached->GetValue(n);
    // first check: only update non-beached particles.
    if (beachedFor < this->beachedAtNumberOfTimes) {
      this->points->GetPoint(n, point);
      // second check: only update points within our grid's boundary.
      if (point[0] <= this->uvGrid->lonMin() || point[0] >= this->uvGrid->lonMax() ||
          point[1] <= this->uvGrid->latMin() || point[1] >= this->uvGrid->latMax()) {
        // sets any particle out of bounds to be beached - so it gets assigned the right colour in the lookup table.
        this->particlesBeached->SetValue(n, this->beachedAtNumberOfTimes);
        continue;
      }

      oldX = point[0];
      oldY = point[1];

      // supersampling
      for (int i = 0; i < SUPERSAMPLINGRATE; i++) {
        int dt = (t - this->lastT) / SUPERSAMPLINGRATE;
        if (dt < 0) {
          // TODO: This is a hack for when the t wraps around,
          // there is probably a more elegant way of dealing with this whole thing
          // that involves having two separate DTs.
          // One for the "render" step time, and one for the computation step time.
          dt = t;
        }
        std::tie(point[1], point[0]) = advector->advect(t, point[1], point[0], dt);

        // Track particle if enabled
        if (isTracking && n == trackedParticleIndex && i == SUPERSAMPLINGRATE - 1) {
            trackedPositions.push_back({point[1], point[0]});
            
            // Calculate velocity *at the new position* after advection
            auto vel = bilinearinterpolate(*uvGrid, t, point[1], point[0]);
            trackedVelocities.push_back({vel.u, vel.v});

            // Use GEBCO-based distance to shore *at the new position* after advection
            double shoreDist = uvGrid->getShoreDistance(point[1], point[0]);
            trackedDistancesToShore.push_back(shoreDist);
        }
        // Track all particles if enabled
        if (isTrackingAll && i == SUPERSAMPLINGRATE - 1) {
            // Add the position *after* advection as the next step's data
            allParticlePositions[n].push_back({point[1], point[0]});
            
            // Calculate velocity *at the new position* after advection
            auto vel = bilinearinterpolate(*uvGrid, t, point[1], point[0]);
            allParticleVelocities[n].push_back({vel.u, vel.v});

            // Use GEBCO-based distance to shore *at the new position* after advection
            double shoreDist = uvGrid->getShoreDistance(point[1], point[0]);
            allParticleDistancesToShore[n].push_back(shoreDist);
        }
      }

      if (this->boundaryType == BoundaryType::Snap) {
          if (isNearestNeighbourZero(*uvGrid, t, point[1], point[0])) {
              this->particlesBeached->SetValue(n, beachedFor + 1);
          } else {
              this->particlesBeached->SetValue(n, std::max(beachedFor - 1, 0));
              this->points->SetPoint(n, point);
              modifiedData = true;
          }
      } else {
          // Logic for FreeSlip and PartialSlip boundary conditions
          // Check if outside grid (permanent beaching)
          if (point[0] <= uvGrid->lonMin() || point[0] >= uvGrid->lonMax() ||
              point[1] <= uvGrid->latMin() || point[1] >= uvGrid->latMax()) {
              this->particlesBeached->SetValue(n, this->beachedAtNumberOfTimes);
          } else {
              // Check distance to shore for beaching criteria
              const double beachingThreshold = 1e-3; // Define a small threshold for "zero" distance
              // Calculate distance to shore after advection for beaching checks
              double shoreDist = uvGrid->getShoreDistance(point[1], point[0]);

              if (shoreDist < beachingThreshold) {
                  // Particle is very close to shore, increment beached count
                  this->particlesBeached->SetValue(n, beachedFor + 1);
                  this->points->SetPoint(n, point); // Update point to new location near shore
                  modifiedData = true;
              } else {
                  // Particle is not beached by distance and within bounds, reset beached status and update point
                  this->particlesBeached->SetValue(n, std::max(beachedFor - 1, 0)); // Decrement, but not below 0
                  this->points->SetPoint(n, point);
                  modifiedData = true;
              }
          }
      }
    }
  }
  if (modifiedData) this->points->Modified();
  this->lastT = t;
}

void LagrangeGlyphs::addObservers(vtkSmartPointer<vtkRenderWindowInteractor> interactor) {
  auto newPointCallBack = createSpawnPointCallback();
  interactor->AddObserver(vtkCommand::LeftButtonPressEvent, newPointCallBack);
  interactor->AddObserver(vtkCommand::LeftButtonReleaseEvent, newPointCallBack);
  interactor->AddObserver(vtkCommand::MouseMoveEvent, newPointCallBack);
  interactor->AddObserver(vtkCommand::RightButtonPressEvent, newPointCallBack);
}

vtkSmartPointer<vtkPoints> LagrangeGlyphs::getPoints() {
  return points;
}

vtkSmartPointer<vtkIntArray> LagrangeGlyphs::getBeached() {
  return particlesBeached;
}

void LagrangeGlyphs::setColour(int red, int green, int blue) {
  actor->GetProperty()->SetColor(red / 255.0, green / 255.0, blue / 255.0);
}

void LagrangeGlyphs::setToDiamond() {
  circleSource->SetGlyphTypeToDiamond();
}

void LagrangeGlyphs::handleGameOver() {
  points->Reset();
  particlesBeached->Reset();
}

void LagrangeGlyphs::startTracking(size_t particleIndex) {
    isTracking = true;
    trackedParticleIndex = particleIndex;
    trackedPositions.clear();
    trackedVelocities.clear();
    trackedDistancesToShore.clear();

    // Add the initial spawn position as the first entry (Step 0)
    if (particleIndex < initialSpawnPositions.size()) {
        trackedPositions.push_back(initialSpawnPositions[particleIndex]);
        // Add dummy velocity and distance to shore for step 0
        trackedVelocities.push_back({0.0, 0.0}); // Assuming zero initial velocity
        double initialShoreDist = uvGrid->getShoreDistance(initialSpawnPositions[particleIndex].first, initialSpawnPositions[particleIndex].second);
        trackedDistancesToShore.push_back(initialShoreDist);
    }
}

void LagrangeGlyphs::stopTracking() {
    isTracking = false;
}

void LagrangeGlyphs::startTrackingAll() {
    isTrackingAll = true;
    isTracking = false;  // Disable single particle tracking
    allParticlePositions.clear();
    allParticleVelocities.clear();
    allParticleDistancesToShore.clear();
    
    // Initialize vectors for each particle with the initial spawn positions at step 0
    size_t numParticles = initialSpawnPositions.size(); // Use size of initialSpawnPositions
    if (numParticles == 0) return; // Nothing to track if no particles spawned

    allParticlePositions.resize(numParticles);
    allParticleVelocities.resize(numParticles);
    allParticleDistancesToShore.resize(numParticles);
    
    for (size_t i = 0; i < numParticles; ++i) {
        allParticlePositions[i].push_back(initialSpawnPositions[i]);
        // Add dummy velocity and distance to shore for step 0
        allParticleVelocities[i].push_back({0.0, 0.0}); // Assuming zero initial velocity
        double initialShoreDist = uvGrid->getShoreDistance(initialSpawnPositions[i].first, initialSpawnPositions[i].second);
        allParticleDistancesToShore[i].push_back(initialShoreDist);
    }
}

void LagrangeGlyphs::stopTrackingAll() {
    isTrackingAll = false;
}

void LagrangeGlyphs::printTrackedParticleInfo(const std::string& outputFilename) const {
    if (!isTracking || trackedPositions.empty()) {
        std::cout << "No particle is being tracked or no data available." << std::endl;
        return;
    }

    // Save to file
    std::string outputPath = "C:/Users/wesle/Documents/Universiteit/2024_2025/Thesis/Opdracht/interactive-track-and-trace/particle-track-and-trace/" + outputFilename;
    std::ofstream outFile(outputPath);
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file for writing: " << outputPath << std::endl;
        return;
    }

    try {
        // Write header
        outFile << "Step,Latitude,Longitude,VelocityU,VelocityV,DistanceToShore\n";
        
        // Write data
        for (size_t i = 0; i < trackedPositions.size(); ++i) {
            const auto& pos = trackedPositions[i];
            const auto& vel = trackedVelocities[i];
            outFile << i << ","
                    << pos.first << ","
                    << pos.second << ","
                    << vel.first << ","
                    << vel.second << ","
                    << trackedDistancesToShore[i] << "\n";
        }
        
        outFile.close();
        std::cout << "\nTrajectory data saved to: " << outputPath << std::endl;
        std::cout << "Run 'python plot_trajectory.py' to visualize the particle path on a map." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Error writing to file: " << e.what() << std::endl;
        if (outFile.is_open()) {
            outFile.close();
        }
    }
}

void LagrangeGlyphs::printAllParticlesInfo(const std::string& outputFilename) const {
    if (!isTrackingAll || allParticlePositions.empty()) {
        std::cout << "No particles are being tracked or no data available." << std::endl;
        return;
    }

    // Save to file
    std::string outputPath = "C:/Users/wesle/Documents/Universiteit/2024_2025/Thesis/Opdracht/interactive-track-and-trace/particle-track-and-trace/" + outputFilename;
    std::ofstream outFile(outputPath);
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file for writing: " << outputPath << std::endl;
        return;
    }

    try {
        // Write header
        outFile << "ParticleID,Step,Latitude,Longitude,VelocityU,VelocityV,DistanceToShore\n";
        
        // Write data for each particle, starting with the initial spawn position at Step 0
        for (size_t particleId = 0; particleId < allParticlePositions.size(); ++particleId) {
            const auto& positions = allParticlePositions[particleId];
            const auto& velocities = allParticleVelocities[particleId];
            const auto& distances = allParticleDistancesToShore[particleId];
            
            for (size_t step = 0; step < positions.size(); ++step) {
                const auto& pos = positions[step];
                const auto& vel = velocities[step];
                outFile << particleId << ","
                        << step << "," // Step index starts from 0 now
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
    } catch (const std::exception& e) {
        std::cerr << "Error writing to file: " << e.what() << std::endl;
        if (outFile.is_open()) {
            outFile.close();
        }
    }
}