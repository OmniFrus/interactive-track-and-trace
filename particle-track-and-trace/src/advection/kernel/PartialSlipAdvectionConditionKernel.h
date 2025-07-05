#ifndef PARTIALSLIPADVECTIONCONDITIONKERNEL_H
#define PARTIALSLIPADVECTIONCONDITIONKERNEL_H

#include "AdvectionKernel.h"
#include "../UVGrid.h"
#include <memory>

class PartialSlipAdvectionConditionKernel : public AdvectionKernel {
public:
    PartialSlipAdvectionConditionKernel(std::unique_ptr<AdvectionKernel> baseKernel, std::shared_ptr<UVGrid> grid);

    std::pair<double, double> advect(int time, double latitude, double longitude, int dt) const override;

private:
    std::unique_ptr<AdvectionKernel> baseKernel;
    std::shared_ptr<UVGrid> grid;
    double slipFactor; // Factor controlling the amount of slip (0 = no slip, 1 = free slip)
};

#endif