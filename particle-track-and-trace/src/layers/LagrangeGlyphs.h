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

/**
 * Implements the Layer class for the case of a Lagrangian visualization.
 * Specifically, this class models the Lagrangian particles in the simulation
 * using the 'glyph' mark and 'transparency' channel to denote age.
 */
class LagrangeGlyphs : public Layer {
public:
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

  /**
   * Load spawn locations from CSV file and create particles
   * @param filename Path to the CSV file containing lat,lon coordinates
   */
  void spawnParticlesFromFile(const std::string& filename);

  vtkSmartPointer<SpawnPointCallback> createSpawnPointCallback();
};

#endif
