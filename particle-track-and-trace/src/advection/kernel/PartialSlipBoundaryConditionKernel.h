#ifndef PARTIALSLIPBOUNDARYCONDITIONKERNEL_H
#define PARTIALSLIPBOUNDARYCONDITIONKERNEL_H

#include "AdvectionKernel.h"
#include "../UVGrid.h"
#include <memory>

class PartialSlipBoundaryConditionKernel : public AdvectionKernel {
public:
    PartialSlipBoundaryConditionKernel(std::unique_ptr<AdvectionKernel> baseKernel, std::shared_ptr<UVGrid> grid);

    std::pair<double, double> advect(int time, double latitude, double longitude, int dt) const override;

private:
    std::unique_ptr<AdvectionKernel> baseKernel;
    std::shared_ptr<UVGrid> grid;
    double slipFactor; // Factor controlling the amount of slip (0 = no slip, 1 = free slip)
};

#endif