#include "FreeSlipBoundaryConditionKernel.h"
#include "../interpolate.h"
#include <cmath>

FreeSlipBoundaryConditionKernel::FreeSlipBoundaryConditionKernel(
    std::unique_ptr<AdvectionKernel> baseKernel,
    std::shared_ptr<UVGrid> grid)
    : baseKernel(std::move(baseKernel))
    , grid(std::move(grid)) {}

std::pair<double, double> FreeSlipBoundaryConditionKernel::advect(int time, double latitude, double longitude, int dt) const {
    // First check if the next position would be on land
    auto [newLat, newLon] = baseKernel->advect(time, latitude, longitude, dt);
    const double epsilon = 1e-5; // About ~1 meter inward (in degrees)

    // Get velocity at current position using bilinear interpolation
    auto vel = bilinearinterpolate(*grid, time, latitude, longitude);
    bool needsReflection = false;

    // Check if the next position would be on land or near coast
    if (isNearestNeighbourZero(*grid, time, newLat, newLon)) {
        needsReflection = true;

        // Calculate normalized coordinates within cell (xi and eta)
        double xi = (longitude - grid->lons[0]) / grid->lonStep();
        xi = xi - std::floor(xi);  // Normalize to [0,1]
        double eta = (latitude - grid->lats[0]) / grid->latStep();
        eta = eta - std::floor(eta);  // Normalize to [0,1]

        const double min_coord = 0.01; // Prevent division by zero

        // Determine which boundary we're closest to and apply appropriate slip
        double distToWest = std::abs(newLon - grid->lonMin());
        double distToEast = std::abs(newLon - grid->lonMax());
        double distToSouth = std::abs(newLat - grid->latMin());
        double distToNorth = std::abs(newLat - grid->latMax());

        // Find closest boundary and apply slip
        if (distToWest <= distToEast && distToWest <= distToSouth && distToWest <= distToNorth) {
            // Western boundary - slip northward/southward
            vel.u = 0;  // Zero normal velocity
            vel.v *= 1.0 / std::max(eta, min_coord);
        } else if (distToEast <= distToSouth && distToEast <= distToNorth) {
            // Eastern boundary - slip northward/southward
            vel.u = 0;  // Zero normal velocity
            vel.v *= 1.0 / std::max(1.0 - eta, min_coord);
        } else if (distToSouth <= distToNorth) {
            // Southern boundary - slip eastward/westward
            vel.v = 0;  // Zero normal velocity
            vel.u *= 1.0 / std::max(xi, min_coord);
        } else {
            // Northern boundary - slip eastward/westward
            vel.v = 0;  // Zero normal velocity
            vel.u *= 1.0 / std::max(1.0 - xi, min_coord);
        }

        // Calculate new position with slip velocities
        newLon = longitude + metreToDegrees(vel.u * dt);
        newLat = latitude + metreToDegrees(vel.v * dt);

        // Ensure we stay within bounds and away from land
        newLat = std::max(grid->latMin() + epsilon, std::min(grid->latMax() - epsilon, newLat));
        newLon = std::max(grid->lonMin() + epsilon, std::min(grid->lonMax() - epsilon, newLon));

        // If new position would still be on land, move back to previous position
        if (isNearestNeighbourZero(*grid, time, newLat, newLon)) {
            newLat = latitude;
            newLon = longitude;
        }
    }

    return {newLat, newLon};
}