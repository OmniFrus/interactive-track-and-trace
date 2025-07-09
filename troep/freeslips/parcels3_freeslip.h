#ifndef PARCELSRK4FREESLIPKERNEL_H
#define PARCELSRK4FREESLIPKERNEL_H

#include "AdvectionKernel.h"
#include "../UVGrid.h"
#include <memory>

class ParcelsRK4FreeSlipKernel : public AdvectionKernel {
public:
    explicit ParcelsRK4FreeSlipKernel(std::shared_ptr<UVGrid> grid);

    std::pair<double, double> advect(int time, double lat, double lon, int dt) const override;

private:
    std::shared_ptr<UVGrid> grid;
};

#endif // PARCELSRK4FREESLIPKERNEL_H