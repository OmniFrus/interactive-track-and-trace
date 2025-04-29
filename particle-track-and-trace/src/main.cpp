#include "layers/BackgroundImage.h"
#include "layers/EulerGlyphs.h"
#include "layers/LagrangeGlyphs.h"
#include "layers/Camera.h"
#include "layers/DayCounter.h"
#include "layers/Timer.h"
#include "Program.h"
#include "advection/UVGrid.h"
#include "advection/kernel/RK4AdvectionKernel.h"
#include "advection/kernel/SnapBoundaryConditionKernel.h"
#include "advection/kernel/FreeSlipBoundaryConditionKernel.h"
#include "advection/kernel/PartialSlipBoundaryConditionKernel.h"
#include "advection/kernel/SnapBoundaryConditionKernel.h"

#include <vtkPolyDataMapper2D.h>
#include <vtkProperty2D.h>
#include <memory>

using namespace std;

constexpr int dt = 60 * 60; // 60 sec/min * 60 mins

int main() {
  cout << "Reading data..." << endl;
  string dataPath = "../../../../../data";
  cout << "Creating UVGrid..." << endl;
  shared_ptr<UVGrid> uvGrid = make_shared<UVGrid>(dataPath);
  cout << "Created UVGrid." << endl;
  
  auto kernelRK4 = make_unique<RK4AdvectionKernel>(uvGrid);
  unique_ptr<AdvectionKernel> boundaryKernel;

    // Choose boundary handling here:
    // For Snap:
    // boundaryKernel = make_unique<SnapBoundaryConditionKernel>(std::move(kernelRK4), uvGrid);
    // Or for Partial Slip:
    // boundaryKernel = make_unique<PartialSlipBoundaryConditionKernel>(std::move(kernelRK4), uvGrid);
    // Or for Free Slip:
     boundaryKernel = make_unique<FreeSlipBoundaryConditionKernel>(std::move(kernelRK4), uvGrid);

  cout << "Starting vtk..." << endl;
  auto program = make_shared<Program>(dt);
  auto timer = make_shared<Timer>(program, dt);
  auto camera = make_shared<Camera>();

  // Create and configure litter particles with spawn locations from CSV
  auto litter = make_shared<LagrangeGlyphs>(uvGrid, std::move(boundaryKernel), dataPath + "/spawn_locations.csv");
  litter->setToDiamond();

  // Create Euler glyphs for flow visualization
  auto euler = make_shared<EulerGlyphs>(uvGrid);
  
  // Create day counter
  auto dayCounter = make_shared<DayCounter>();

  // Add layers to program in correct order
  program->addLayer(timer);
  program->addLayer(make_shared<BackgroundImage>(dataPath + "/northsea.png"));
  program->addLayer(euler);  // Add Euler glyphs before litter
  program->addLayer(litter);
  program->addLayer(camera);
  program->addLayer(dayCounter);

  program->render();

  return EXIT_SUCCESS;
}

