#include "PartialSlipBoundaryConditionKernel.h"
#include "../interpolate.h"
#include <cmath>
#include <algorithm>

PartialSlipBoundaryConditionKernel::PartialSlipBoundaryConditionKernel(
    std::unique_ptr<AdvectionKernel> baseKernel,
    std::shared_ptr<UVGrid> grid)
    : baseKernel(std::move(baseKernel))
    , grid(std::move(grid)) {}

std::pair<double, double> PartialSlipBoundaryConditionKernel::advect(int time, double latitude, double longitude, int dt) const {
    auto [newLat, newLon] = baseKernel->advect(time, latitude, longitude, dt);

    const double epsilon = 1e-5;
    const double shoreThreshold = 5000.0; // 5 km
    const double minVelocity = 1e-4;
    const double maxVelocity = 3.0;
    const double minCoord = 0.05;
    const double slipRatio = 0.5;  // Between 0 (no slip) and 1 (full slip)

    auto vel = bilinearinterpolate(*grid, time, latitude, longitude);

    if (grid->isNearShore(latitude, longitude, shoreThreshold)) {
        double xi = (longitude - grid->lonMin()) / grid->lonStep();
        xi = std::clamp(xi - std::floor(xi), minCoord, 1.0 - minCoord);
        double eta = (latitude - grid->latMin()) / grid->latStep();
        eta = std::clamp(eta - std::floor(eta), minCoord, 1.0 - minCoord);

        double distToWest  = std::abs(newLon - grid->lonMin());
        double distToEast  = std::abs(newLon - grid->lonMax());
        double distToSouth = std::abs(newLat - grid->latMin());
        double distToNorth = std::abs(newLat - grid->latMax());

        if (distToWest <= distToEast && distToWest <= distToSouth && distToWest <= distToNorth) {
            vel.v *= (1.0 - 1.0 / (2.0 * xi)) / (1.0 - xi);
        } else if (distToEast <= distToSouth && distToEast <= distToNorth) {
            vel.v *= (1.0 - 1.0 / (2.0 * (1.0 - xi)) / xi);
        } else if (distToSouth <= distToNorth) {
            vel.u *= (1.0 - 1.0 / (2.0 * eta)) / (1.0 - eta);
        } else {
            vel.u *= (1.0 - 1.0 / (2.0 * (1.0 - eta)) / eta);
        }

        vel.u = std::clamp(vel.u, -maxVelocity, maxVelocity);
        vel.v = std::clamp(vel.v, -maxVelocity, maxVelocity);

        newLon = longitude + metreToDegrees(vel.u * dt);
        newLat = latitude + metreToDegrees(vel.v * dt);

        newLat = std::clamp(newLat, grid->latMin() + epsilon, grid->latMax() - epsilon);
        newLon = std::clamp(newLon, grid->lonMin() + epsilon, grid->lonMax() - epsilon);
    }

    return {newLat, newLon};
}
