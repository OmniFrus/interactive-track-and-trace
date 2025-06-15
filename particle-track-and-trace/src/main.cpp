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
#include "advection/kernel/ParcelsBoundaryConditionKernel.h"

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

  // Choose boundary condition here:
  // For Snap:
  // boundaryKernel = make_unique<SnapBoundaryConditionKernel>(std::move(kernelRK4), uvGrid);
  // For Partial Slip:
  // boundaryKernel = make_unique<PartialSlipBoundaryConditionKernel>(std::move(kernelRK4), uvGrid);
  // For Free Slip:
  // boundaryKernel = make_unique<FreeSlipBoundaryConditionKernel>(std::move(kernelRK4), uvGrid);
  // For Parcels Slip:
   boundaryKernel = make_unique<ParcelsRK4FreeSlipKernel>(std::move(kernelRK4),uvGrid);
  cout << "Created boundaryKernel successfully." << endl;

  cout << "Starting vtk..." << endl;
  auto program = make_shared<Program>(dt);
  cout << "Created Program." << endl;
  
  auto timer = make_shared<Timer>(program, dt);
  cout << "Created Timer." << endl;
  
  auto camera = make_shared<Camera>();
  cout << "Created Camera." << endl;

  // Create and configure litter particles with spawn locations from CSV
  auto litter = make_shared<LagrangeGlyphs>(uvGrid, std::move(boundaryKernel), dataPath + "/spawn_locations.csv");
  cout << "Created LagrangeGlyphs successfully." << endl;
  
  litter->setToDiamond();
  litter->setCoastalTimeThreshold(6); // Set to 6 hours
  
  // Choose beaching conditions here:
  // litter->setBeachingType(LagrangeGlyphs::BeachingType::VelocityBased);    // Original snap boundary logic
  // litter->setBeachingType(LagrangeGlyphs::BeachingType::DistanceBased);    // Based on distance to shore
   litter->setBeachingType(LagrangeGlyphs::BeachingType::DirectionalBased); // Based on direction and distance (dual condition)
  // litter->setBeachingType(LagrangeGlyphs::BeachingType::None);              // No beaching (only out of bounds is considered beached)

  // Enable/disable directional check for DirectionalBased beaching
  litter->setEnableDirectionalCheck(true);

  // If u want to track all particles, use this:
   litter->startTrackingAll();

  // If u only want to track a single particle, use this:
   litter->startTracking(100);

  // Create Euler glyphs for flow visualization
  auto euler = make_shared<EulerGlyphs>(uvGrid);
  cout << "Created EulerGlyphs successfully." << endl;

  // Create day counter
  auto dayCounter = make_shared<DayCounter>();
  cout << "Created DayCounter successfully." << endl;

  // Add layers to program in correct order
  program->addLayer(timer);
  cout << "Added Timer layer." << endl;
  
  program->addLayer(make_shared<BackgroundImage>(dataPath + "/northsea.png"));
  cout << "Added BackgroundImage layer." << endl;

  program->addLayer(euler);  // Add Euler glyphs before litter
  cout << "Added EulerGlyphs layer." << endl;
  
  program->addLayer(litter);
  cout << "Added LagrangeGlyphs layer." << endl;
  
  program->addLayer(camera);
  cout << "Added Camera layer." << endl;
  
  program->addLayer(dayCounter);
  cout << "Added DayCounter layer." << endl;

  cout << "Starting render..." << endl;
  program->render();
  cout << "Render completed." << endl;

  // Automatically print tracking information based on active tracking mode
  if (litter->isTrackingAll()) {
    litter->printAllParticlesInfo("all_particles_trajectory.csv");
  }
  if (litter->isTracking()) {
    litter->printTrackedParticleInfo("single_particle_trajectory.csv");
  }

  return EXIT_SUCCESS;
}

