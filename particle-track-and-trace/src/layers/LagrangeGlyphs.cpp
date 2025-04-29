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
      }

      // if the particle's location remains unchanged, increase beachedFor number. Else, decrease it and update point position.
      //if (abs(oldX - point[0]) < EPS and abs(oldY-point[1]) < EPS) {

      bool useVelocityBeaching = (this->boundaryType == BoundaryType::Snap);
      if (useVelocityBeaching) {
          if (isNearestNeighbourZero(*uvGrid, t, point[1], point[0])) {
              this->particlesBeached->SetValue(n, beachedFor + 1);
          } else {
              this->particlesBeached->SetValue(n, std::max(beachedFor - 1, 0));
              this->points->SetPoint(n, point);
              modifiedData = true;
          }
      } else {
          // In FreeSlip or PartialSlip: only beach if outside grid
          if (point[0] <= uvGrid->lonMin() || point[0] >= uvGrid->lonMax() ||
              point[1] <= uvGrid->latMin() || point[1] >= uvGrid->latMax()) {
              this->particlesBeached->SetValue(n, this->beachedAtNumberOfTimes);
          } else {
              this->particlesBeached->SetValue(n, 0);
              this->points->SetPoint(n, point);
              modifiedData = true;
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
