#ifndef PARCELSSLIPADVECTIONKERNEL_H
#define PARCELSSLIPADVECTIONKERNEL_H

#include "AdvectionKernel.h"
#include "RK4AdvectionKernel.h"
#include "../UVGrid.h"
#include <memory>

class ParcelsAdvectionKernel : public AdvectionKernel {
public:
    ParcelsAdvectionKernel(std::unique_ptr<RK4AdvectionKernel> baseKernel, std::shared_ptr<UVGrid> grid);

    std::pair<double, double> advect(int time, double lat, double lon, int dt) const override;

private:
    std::unique_ptr<RK4AdvectionKernel> baseKernel;
    std::shared_ptr<UVGrid> grid;
};

#endif // PARCELSRK4FREESLIPKERNEL_H