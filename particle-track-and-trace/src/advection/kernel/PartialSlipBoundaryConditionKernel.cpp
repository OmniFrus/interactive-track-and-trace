#include "PartialSlipBoundaryConditionKernel.h"

PartialSlipBoundaryConditionKernel::PartialSlipBoundaryConditionKernel(std::unique_ptr<AdvectionKernel> baseKernel, std::shared_ptr<UVGrid> grid)
    : baseKernel(std::move(baseKernel)), grid(std::move(grid)) {}

std::pair<double, double> PartialSlipBoundaryConditionKernel::advect(int time, double latitude, double longitude, int dt) const {
    auto [newLat, newLon] = baseKernel->advect(time, latitude, longitude, dt);

    const double epsilon = 1e-5; // About ~1 meter inward

    // Stop movement across boundary (freeze component)
    if (newLat < grid->latMin() || newLat > grid->latMax()) {
        newLat = latitude;
    }
    if (newLon < grid->lonMin() || newLon > grid->lonMax()) {
        newLon = longitude;
    }

    // Push slightly inward to avoid sticking exactly at the border
    if (newLat < grid->latMin()) newLat = grid->latMin() + epsilon;
    if (newLat > grid->latMax()) newLat = grid->latMax() - epsilon;
    if (newLon < grid->lonMin()) newLon = grid->lonMin() + epsilon;
    if (newLon > grid->lonMax()) newLon = grid->lonMax() - epsilon;

    return {newLat, newLon};
}