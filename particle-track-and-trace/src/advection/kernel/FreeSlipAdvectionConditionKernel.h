#ifndef FREESLIPADVECTIONCONDITIONKERNEL_H
#define FREESLIPADVECTIONCONDITIONKERNEL_H

#include "AdvectionKernel.h"
#include "../UVGrid.h"
#include <memory>

class FreeSlipAdvectionConditionKernel : public AdvectionKernel {
public:
    FreeSlipAdvectionConditionKernel(std::unique_ptr<AdvectionKernel> baseKernel, std::shared_ptr<UVGrid> grid);

    std::pair<double, double> advect(int time, double latitude, double longitude, int dt) const override;

private:
    std::unique_ptr<AdvectionKernel> baseKernel;
    std::shared_ptr<UVGrid> grid;
};

#endif