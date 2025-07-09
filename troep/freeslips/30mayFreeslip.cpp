#include "FreeSlipBoundaryConditionKernel.h"
#include "../interpolate.h"
#include <cmath>
#include <algorithm>

FreeSlipBoundaryConditionKernel::FreeSlipBoundaryConditionKernel(
    std::unique_ptr<AdvectionKernel> baseKernel,
    std::shared_ptr<UVGrid> grid)
    : baseKernel(std::move(baseKernel))
    , grid(std::move(grid)) {}

std::pair<double, double> FreeSlipBoundaryConditionKernel::advect(int time, double latitude, double longitude, int dt) const {
    auto [newLat, newLon] = baseKernel->advect(time, latitude, longitude, dt);
    const double epsilon = 1e-5;
    const double shoreThreshold = 2000.0; // 2km
    const double minVelocity = 1e-4;
    const double maxVelocity = 3.0;
    const double minCoord = 0.05;
    const double displacementStrength = 1.0; // Strength of the displacement field

    auto vel = bilinearinterpolate(*grid, time, latitude, longitude);

    // Get current shore distance
    double shoreDist = grid->getShoreDistance(latitude, longitude);
    
    // Calculate displacement field effect
    if (shoreDist < shoreThreshold) {
        // Calculate normalized coordinates within cell
        double xi = (longitude - grid->lonMin()) / grid->lonStep();
        xi = std::clamp(xi - std::floor(xi), minCoord, 1.0 - minCoord);
        double eta = (latitude - grid->latMin()) / grid->latStep();
        eta = std::clamp(eta - std::floor(eta), minCoord, 1.0 - minCoord);

        // Calculate distance to boundaries
        double distToWest  = std::abs(newLon - grid->lonMin());
        double distToEast  = std::abs(newLon - grid->lonMax());
        double distToSouth = std::abs(newLat - grid->latMin());
        double distToNorth = std::abs(newLat - grid->latMax());

        // Apply free-slip boundary condition
        if (distToWest <= distToEast && distToWest <= distToSouth && distToWest <= distToNorth) {
            vel.u = 0;
            vel.v *= 1.0 / xi; // Formula 4
        } else if (distToEast <= distToSouth && distToEast <= distToNorth) {
            // Closest to East (vertical boundary): normal vel.u = 0, apply xi-formula to vel.v
            vel.u = 0;
            vel.v *= 1.0 / (1.0 - xi); // Formula 8
        } else if (distToSouth <= distToNorth) {
            // Closest to South (horizontal boundary): normal vel.v = 0, apply eta-formula to vel.u
            vel.v = 0;
            vel.u *= 1.0 / eta; // Formula 1
        } else {
            // Closest to North (horizontal boundary): normal vel.v = 0, apply eta-formula to vel.u
            vel.v = 0;
            vel.u *= 1.0 / (1.0 - eta); // Formula 2
        }

        // Add displacement field effect
        double displacementFactor = (1.0 - shoreDist / shoreThreshold) * displacementStrength;
        
        // Calculate displacement direction (away from shore)
        double dx = 0.0, dy = 0.0;
        if (distToWest <= distToEast && distToWest <= distToSouth && distToWest <= distToNorth) {
            dx = displacementFactor; // Move east
        } else if (distToEast <= distToSouth && distToEast <= distToNorth) {
            dx = -displacementFactor; // Move west
        } else if (distToSouth <= distToNorth) {
            dy = displacementFactor; // Move north
        } else {
            dy = -displacementFactor; // Move south
        }

        // Apply displacement
        vel.u += dx;
        vel.v += dy;

        // Clamp velocities
        vel.u = std::clamp(vel.u, -maxVelocity, maxVelocity);
        vel.v = std::clamp(vel.v, -maxVelocity, maxVelocity);

        // Update position
        newLon = longitude + metreToDegrees(vel.u * dt);
        newLat = latitude + metreToDegrees(vel.v * dt);

        // Ensure we stay within bounds
        newLat = std::clamp(newLat, grid->latMin() + epsilon, grid->latMax() - epsilon);
        newLon = std::clamp(newLon, grid->lonMin() + epsilon, grid->lonMax() - epsilon);
    }

    return {newLat, newLon};
}
