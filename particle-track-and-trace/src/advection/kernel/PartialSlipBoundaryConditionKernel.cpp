#include "PartialSlipBoundaryConditionKernel.h"
#include "../interpolate.h"
#include <cmath>

PartialSlipBoundaryConditionKernel::PartialSlipBoundaryConditionKernel(
    std::unique_ptr<AdvectionKernel> baseKernel, 
    std::shared_ptr<UVGrid> grid)
    : baseKernel(std::move(baseKernel))
    , grid(std::move(grid))
    , slipFactor(0.5) {} // Default to 0.5 for partial slip

std::pair<double, double> PartialSlipBoundaryConditionKernel::advect(int time, double latitude, double longitude, int dt) const {
    // First check if the next position would be on land
    auto [newLat, newLon] = baseKernel->advect(time, latitude, longitude, dt);
    const double epsilon = 1e-5; // About ~1 meter inward (in degrees)
    const double shoreThreshold = 1000.0; // 1km threshold for near-shore behavior

    // Get velocity at current position using bilinear interpolation
    auto vel = bilinearinterpolate(*grid, time, latitude, longitude);
    bool needsReflection = false;

    // Find grid indices for current position
    size_t latIndex = static_cast<size_t>((latitude - grid->latMin()) / grid->latStep());
    size_t lonIndex = static_cast<size_t>((longitude - grid->lonMin()) / grid->lonStep());

    // Check if we're near shore or if next position would be on land
    if (grid->isNearShore(latIndex, lonIndex, shoreThreshold) || 
        isNearestNeighbourZero(*grid, time, newLat, newLon)) {
        needsReflection = true;

        // Calculate normalized coordinates within cell (xi and eta)
        double xi = (longitude - grid->lons[0]) / grid->lonStep();
        xi = xi - std::floor(xi);  // Normalize to [0,1]
        double eta = (latitude - grid->lats[0]) / grid->latStep();
        eta = eta - std::floor(eta);  // Normalize to [0,1]

        const double min_coord = 0.01; // Prevent division by zero

        // Get distance to shore at current position
        double distToShore = grid->getShoreDistance(latIndex, lonIndex);
        
        // Scale slip factor based on distance to shore
        double shoreSlipFactor = std::min(1.0, distToShore / shoreThreshold);
        double effectiveSlipFactor = slipFactor * shoreSlipFactor;

        // Determine which boundary we're closest to and apply appropriate slip
        double distToWest = std::abs(newLon - grid->lonMin());
        double distToEast = std::abs(newLon - grid->lonMax());
        double distToSouth = std::abs(newLat - grid->latMin());
        double distToNorth = std::abs(newLat - grid->latMax());

        // Find closest boundary and apply partial slip factors
        if (distToWest <= distToEast && distToWest <= distToSouth && distToWest <= distToNorth) {
            // Western boundary - partial slip northward/southward
            vel.u = 0;  // Zero normal velocity
            // f_v = (1/2 + 1/2η)/η for western boundary
            vel.v *= effectiveSlipFactor * (0.5 + 0.5 * eta) / std::max(eta, min_coord);
        } else if (distToEast <= distToSouth && distToEast <= distToNorth) {
            // Eastern boundary - partial slip northward/southward
            vel.u = 0;  // Zero normal velocity
            // f_v = (1 - 1/2η)/(1-η) for eastern boundary
            vel.v *= effectiveSlipFactor * (1.0 - 0.5 * eta) / std::max(1.0 - eta, min_coord);
        } else if (distToSouth <= distToNorth) {
            // Southern boundary - partial slip eastward/westward
            vel.v = 0;  // Zero normal velocity
            // f_u = (1/2 + 1/2ξ)/ξ for southern boundary
            vel.u *= effectiveSlipFactor * (0.5 + 0.5 * xi) / std::max(xi, min_coord);
        } else {
            // Northern boundary - partial slip eastward/westward
            vel.v = 0;  // Zero normal velocity
            // f_u = (1 - 1/2ξ)/(1-ξ) for northern boundary
            vel.u *= effectiveSlipFactor * (1.0 - 0.5 * xi) / std::max(1.0 - xi, min_coord);
        }

        // Calculate new position with partial slip velocities
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