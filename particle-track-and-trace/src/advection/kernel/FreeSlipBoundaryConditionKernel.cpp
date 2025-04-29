#include "FreeSlipBoundaryConditionKernel.h"

FreeSlipBoundaryConditionKernel::FreeSlipBoundaryConditionKernel(std::unique_ptr<AdvectionKernel> baseKernel, std::shared_ptr<UVGrid> grid)
    : baseKernel(std::move(baseKernel)), grid(std::move(grid)) {}

std::pair<double, double> FreeSlipBoundaryConditionKernel::advect(int time, double latitude, double longitude, int dt) const {
    auto [newLat, newLon] = baseKernel->advect(time, latitude, longitude, dt);

    const double epsilon = 1e-5; // About ~1 meter inward (in degrees)

    // Reflect if out of bounds
    if (newLat < grid->latMin()) newLat = 2 * grid->latMin() - newLat;
    if (newLat > grid->latMax()) newLat = 2 * grid->latMax() - newLat;
    if (newLon < grid->lonMin()) newLon = 2 * grid->lonMin() - newLon;
    if (newLon > grid->lonMax()) newLon = 2 * grid->lonMax() - newLon;

    // After reflection, push slightly inward to avoid boundary sticking
    if (newLat < grid->latMin()) newLat = grid->latMin() + epsilon;
    if (newLat > grid->latMax()) newLat = grid->latMax() - epsilon;
    if (newLon < grid->lonMin()) newLon = grid->lonMin() + epsilon;
    if (newLon > grid->lonMax()) newLon = grid->lonMax() - epsilon;

    return {newLat, newLon};
}