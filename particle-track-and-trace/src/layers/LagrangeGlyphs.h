#ifndef LAGRANGEGLYPHS_H
#define LAGRANGEGLYPHS_H

#include "Layer.h"
#include "../advection/kernel/AdvectionKernel.h"
#include "../commands/SpawnPointCallback.h"
#include "../gameovers/GameoverCallback.h"
#include <vtkPolyData.h>
#include <vtkInteractorStyle.h>
#include <vtkGlyphSource2D.h>
#include <string>
#include <vector>
#include <utility>

/**
 * Implements the Layer class for the case of a Lagrangian visualization.
 * Specifically, this class models the Lagrangian particles in the simulation
 * using the 'glyph' mark and 'transparency' channel to denote age.
 */
class LagrangeGlyphs : public Layer {
public:
  enum class BeachingType {
    VelocityBased,    // Original snap boundary logic
    DistanceBased,    // Based on distance to shore
    DirectionalBased, // Based on direction and distance
    None             // No beaching (like FreeSlip)
  };

  /**
   * Constructor.
   * @param uvGrid UVGrid used for boundary conditions calculations
   * @param advectionKernel advects particles using given kernel
   * @param spawnFile Path to CSV file containing spawn locations
   */
  LagrangeGlyphs(std::shared_ptr<UVGrid> uvGrid, std::unique_ptr<AdvectionKernel> advectionKernel, const std::string& spawnFile);
  /**
   * This function spoofs a few points in the dataset. Mostly used for testing.
   */
  void spoofPoints();

  /**
   * updates the glyphs to reflect the given timestamp in the dataset.
   * @param t : the time at which to fetch the data.
   */
  void updateData(int t) override;

  void addObservers(vtkSmartPointer<vtkRenderWindowInteractor> interactor) override;

  vtkSmartPointer<vtkPoints> getPoints();

  vtkSmartPointer<vtkIntArray> getBeached();

  void setColour(int red, int green, int blue);

  void setToDiamond();

  void handleGameOver() override;

  // Directional check control
  void setEnableDirectionalCheck(bool enabled);


  // Tracking methods
  void startTracking(size_t particleIndex);
  void stopTracking();
  void startTrackingAll();
  void stopTrackingAll();
  void printTrackedParticleInfo(const std::string& outputFilename) const;
  void printAllParticlesInfo(const std::string& outputFilename) const;
  const std::vector<std::pair<double, double>>& getTrackedPositions() const { return trackedPositions; }
  const std::vector<std::pair<double, double>>& getTrackedVelocities() const { return trackedVelocities; }
  const std::vector<double>& getTrackedDistancesToShore() const { return trackedDistancesToShore; }
  const std::vector<int>& getCoastalResidenceTimes() const { return coastalResidenceTimes; }

  // Set the minimum coastal residence time (in hours) required for beaching
  void setCoastalTimeThreshold(double hours);

  // Add new method to set beaching type
  void setBeachingType(BeachingType type) { beachingType = type; }

private:
  vtkNew<vtkPoints> points;
  vtkSmartPointer<vtkPolyData> data;
  vtkSmartPointer<vtkIntArray> particlesBeached;
  std::unique_ptr<AdvectionKernel> advector;
  std::shared_ptr<UVGrid> uvGrid;
  vtkNew<vtkActor> actor;
  vtkNew<vtkGlyphSource2D> circleSource;
  int lastT = 1000;
  int beachedAtNumberOfTimes = 50;

  // Tracking variables
  bool isTracking = false;
  bool isTrackingAll = false;
  size_t trackedParticleIndex = 0;
  std::vector<std::pair<double, double>> trackedPositions;
  std::vector<std::pair<double, double>> trackedVelocities;
  std::vector<double> trackedDistancesToShore;
  
  // All particles tracking
  std::vector<std::vector<std::pair<double, double>>> allParticlePositions;
  std::vector<std::vector<std::pair<double, double>>> allParticleVelocities;
  std::vector<std::vector<double>> allParticleDistancesToShore;
  std::vector<int> coastalResidenceTimes;
  std::vector<int> dualConditionResidenceTimes;
  int coastalTimeThreshold = 6 * 3600; // minimum seconds in coastal buffer before beaching
  std::vector<std::pair<double, double>> initialSpawnPositions; // Store initial positions
  bool enableDirectionalCheck = true; // Direction check for beaching

  BeachingType beachingType = BeachingType::VelocityBased; // Default to velocity based

  /**
   * Load spawn locations from CSV file and create particles
   * @param filename Path to the CSV file containing lat,lon coordinates
   */
  void spawnParticlesFromFile(const std::string& filename);

  vtkSmartPointer<SpawnPointCallback> createSpawnPointCallback();

  enum class BoundaryType {
      Snap,
      FreeSlip,
      PartialSlip,
      DualCondition
  };
  
  BoundaryType boundaryType;
};

#endif
