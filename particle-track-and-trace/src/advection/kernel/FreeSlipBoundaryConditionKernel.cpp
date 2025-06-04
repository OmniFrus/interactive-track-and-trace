#include "FreeSlipBoundaryConditionKernel.h"
#include "../interpolate.h"
#include <cmath>
#include <algorithm>
#include <iostream>

FreeSlipBoundaryConditionKernel::FreeSlipBoundaryConditionKernel(
    std::unique_ptr<AdvectionKernel> baseKernel,
    std::shared_ptr<UVGrid> grid)
    : baseKernel(std::move(baseKernel))
    , grid(std::move(grid)) {}

std::pair<double, double> FreeSlipBoundaryConditionKernel::advect(int time, double latitude, double longitude, int dt) const {
    auto [newLat, newLon] = baseKernel->advect(time, latitude, longitude, dt);

    const double epsilon = 1e-5;
    const double shoreThreshold = 5000.0; // 2 km
    const double minVelocity = 1e-4;
    const double maxVelocity = 3.0;
    const double minCoord = 0.05;

    // Get velocity at the current position (needed for slip logic)
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

        // Apply Parcels-style free slip based on boundary direction
        // Using their exact scaling factors:
        // 1: f_u = 1/η, η = eta 
        // 2: f_u = 1/(1-η)
        // 4: f_v = 1/ξ, ξ = xi
        // 8: f_v = 1/(1-ξ)
        if (distToWest <= distToEast && distToWest <= distToSouth && distToWest <= distToNorth) {
            // Western boundary - slip northward/southward
            vel.u = 0;  // Zero normal velocity
            vel.v *= 1.0 / xi;  // Formula 4
        } else if (distToEast <= distToSouth && distToEast <= distToNorth) {
            // Eastern boundary - slip northward/southward
            vel.u = 0;  // Zero normal velocity
            vel.v *= 1.0 / (1.0 - xi);  // Formula 8
        } else if (distToSouth <= distToNorth) {
            // Southern boundary - slip eastward/westward
            vel.v = 0;  // Zero normal velocity
            vel.u *= 1.0 / eta;  // Formula 1 
        } else {
            // Northern boundary - slip eastward/westward
            vel.v = 0;  // Zero normal velocity
            vel.u *= 1.0 / (1.0 - eta);  // Formula 2
        }

        vel.u = std::clamp(vel.u, -maxVelocity, maxVelocity);
        vel.v = std::clamp(vel.v, -maxVelocity, maxVelocity);

        newLon = longitude + metreToDegrees(vel.u * dt);
        newLat = latitude + metreToDegrees(vel.v * dt);

        newLat = std::clamp(newLat, grid->latMin() + epsilon, grid->latMax() - epsilon);
        newLon = std::clamp(newLon, grid->lonMin() + epsilon, grid->lonMax() - epsilon);
    }
    // Reject motion into coast
    double newShoreDist = grid->getShoreDistance(newLat, newLon);
    if (newShoreDist < 200.0) {  // Within 100m of land
        return {latitude, longitude};  // Stay in current location
    }

    return {newLat, newLon};
}