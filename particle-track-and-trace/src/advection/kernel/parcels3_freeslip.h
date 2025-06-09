#ifndef PARCELSRK4FREESLIPKERNEL_H
#define PARCELSRK4FREESLIPKERNEL_H

#include "AdvectionKernel.h"
#include "RK4AdvectionKernel.h"
#include "../UVGrid.h"
#include <memory>

class ParcelsRK4FreeSlipKernel : public AdvectionKernel {
public:
    ParcelsRK4FreeSlipKernel(std::unique_ptr<RK4AdvectionKernel> baseKernel, std::shared_ptr<UVGrid> grid);

    std::pair<double, double> advect(int time, double lat, double lon, int dt) const override;

private:
    std::unique_ptr<RK4AdvectionKernel> baseKernel;
    std::shared_ptr<UVGrid> grid;
};

#endif // PARCELSRK4FREESLIPKERNEL_H