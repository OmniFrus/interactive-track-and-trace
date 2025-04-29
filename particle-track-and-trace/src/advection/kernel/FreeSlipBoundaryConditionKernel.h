#ifndef FREESLIPBOUNDARYCONDITIONKERNEL_H
#define FREESLIPBOUNDARYCONDITIONKERNEL_H

#include "AdvectionKernel.h"
#include "../UVGrid.h"
#include <memory>

class FreeSlipBoundaryConditionKernel : public AdvectionKernel {
public:
    FreeSlipBoundaryConditionKernel(std::unique_ptr<AdvectionKernel> baseKernel, std::shared_ptr<UVGrid> grid);

    std::pair<double, double> advect(int time, double latitude, double longitude, int dt) const override;

private:
    std::unique_ptr<AdvectionKernel> baseKernel;
    std::shared_ptr<UVGrid> grid;
};

#endif