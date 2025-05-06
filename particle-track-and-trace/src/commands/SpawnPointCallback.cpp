#include "SpawnPointCallback.h"

#include <vtkVertex.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkCommand.h>
#include <vtkRenderWindow.h>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <algorithm>

#include "../CartographicTransformation.h"

void convertDisplayToWorld(vtkRenderer *renderer, int x, int y, double *worldPos) {
    double displayPos[3] = {static_cast<double>(x), static_cast<double>(y), 0.0};
    renderer->SetDisplayPoint(displayPos);
    renderer->DisplayToWorld();
    renderer->GetWorldPoint(worldPos);
}

void SpawnPointCallback::Execute(vtkObject *caller, unsigned long evId, void *callData) {
    auto interactor = reinterpret_cast<vtkRenderWindowInteractor *>(caller);

    int x, y;
    interactor->GetEventPosition(x, y);

    double displayPos[3] = {static_cast<double>(x), static_cast<double>(y), 0.0};
    double worldPos[4] = {0, 0, 0, 0};
    ren->SetDisplayPoint(displayPos);
    ren->DisplayToWorld();
    ren->GetWorldPoint(worldPos);
    inverseCartographicProjection->TransformPoint(worldPos, worldPos);

    // LEFT CLICK = SPAWN PARTICLE
    if (evId == vtkCommand::LeftButtonPressEvent) {
        std::cout << "spawned at " << worldPos[1] << ", " << worldPos[0] << std::endl;

        points->InsertNextPoint(worldPos[0], worldPos[1], 0);
        particlesBeached->InsertNextValue(0);

        std::ofstream outFile("../../../../../data/spawn_locations.csv", std::ios::app);
        if (outFile.is_open()) {
            outFile << std::fixed << std::setprecision(6)
                    << worldPos[1] << "," << worldPos[0] << "\n";  // lat, lon
            outFile.close();
        } else {
            std::cerr << "Unable to open file for writing." << std::endl;
        }

        this->particlesBeached->Modified();
        this->points->Modified();
        ren->GetRenderWindow()->Render();
        return;
    }

    // RIGHT CLICK = DELETE PARTICLE
    if (evId == vtkCommand::RightButtonPressEvent) {
        for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i) {
            double p[3];
            points->GetPoint(i, p);
            if (std::abs(p[0] - worldPos[0]) < 0.1 && std::abs(p[1] - worldPos[1]) < 0.1) {
                // Rebuild VTK points and beaching array without the deleted point
                vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
                vtkSmartPointer<vtkIntArray> newBeached = vtkSmartPointer<vtkIntArray>::New();
                newBeached->SetNumberOfComponents(1);

                for (vtkIdType j = 0; j < points->GetNumberOfPoints(); ++j) {
                    if (j == i) continue;
                    double temp[3];
                    points->GetPoint(j, temp);
                    newPoints->InsertNextPoint(temp);
                    newBeached->InsertNextValue(particlesBeached->GetValue(j));
                }

                points->DeepCopy(newPoints);
                particlesBeached->DeepCopy(newBeached);

                // CSV update
                std::vector<std::string> lines;
                std::ifstream inFile("../../../../../data/spawn_locations.csv");
                std::string line;
                std::string header;

                // Save header separately
                if (std::getline(inFile, line)) {
                    header = line;
                    lines.push_back(header);
                }

                while (std::getline(inFile, line)) {
                    lines.push_back(line);
                }
                inFile.close();

                // Delete the matching lat,lon line using coordinate comparison with tolerance
                double lat = p[1], lon = p[0];
                auto it = std::remove_if(lines.begin() + 1, lines.end(), [&](const std::string &line) {
                    std::stringstream ss(line);
                    std::string lat_str, lon_str;
                    if (std::getline(ss, lat_str, ',') && std::getline(ss, lon_str, ',')) {
                        try {
                            double csvLat = std::stod(lat_str);
                            double csvLon = std::stod(lon_str);
                            return std::abs(csvLat - lat) < 1e-5 && std::abs(csvLon - lon) < 1e-5;
                        } catch (...) {
                            return false;
                        }
                    }
                    return false;
                });
                lines.erase(it, lines.end());

                // Rewrite CSV file
                std::ofstream outFileUpdate("../../../../../data/spawn_locations.csv");
                for (const auto &l : lines) {
                    outFileUpdate << l << "\n";
                }
                outFileUpdate.close();

                this->particlesBeached->Modified();
                this->points->Modified();
                ren->GetRenderWindow()->Render();
                return;
            }
        }
    }
}

SpawnPointCallback::SpawnPointCallback() : data(nullptr),
                                           points(nullptr),
                                           inverseCartographicProjection(nullptr),
                                           uvGrid(nullptr) {}

SpawnPointCallback *SpawnPointCallback::New() {
    return new SpawnPointCallback;
}

void SpawnPointCallback::setData(const vtkSmartPointer<vtkPolyData> &data) {
    this->data = data;
}

void SpawnPointCallback::setPoints(const vtkSmartPointer<vtkPoints> &points) {
    this->points = points;
}

void SpawnPointCallback::setRen(const vtkSmartPointer<vtkRenderer> &ren) {
    this->ren = ren;
}

void SpawnPointCallback::setUVGrid(const std::shared_ptr<UVGrid> &uvGrid) {
    this->uvGrid = uvGrid;
    inverseCartographicProjection = createInverseCartographicTransformFilter(*uvGrid)->GetTransform();
}

void SpawnPointCallback::setBeached(const vtkSmartPointer<vtkIntArray> &ints) {
    this->particlesBeached = ints;
}
