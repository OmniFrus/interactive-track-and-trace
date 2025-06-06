#ifndef DUALBOUNDARYCONDITIONKERNEL_H
#define DUALBOUNDARYCONDITIONKERNEL_H

#include "AdvectionKernel.h"
#include "../UVGrid.h"
#include <memory>
#include <algorithm>

/**
 * Implements a boundary handling strategy that only beaches particles
 * when they move landward into the coastline. A coastal buffer zone
 * of 2 km around land is used to evaluate the motion relative to the
 * shoreline based on the gradient of the distance-to-shore field.
 */
class DualBoundaryConditionKernel : public AdvectionKernel {
public:
    DualBoundaryConditionKernel(std::unique_ptr<AdvectionKernel> baseKernel, std::shared_ptr<UVGrid> grid);

    std::pair<double, double> advect(int time, double latitude, double longitude, int dt) const override;

private:
    std::unique_ptr<AdvectionKernel> baseKernel;
    std::shared_ptr<UVGrid> grid;
};

#endif //DUALBOUNDARYCONDITIONKERNEL_H