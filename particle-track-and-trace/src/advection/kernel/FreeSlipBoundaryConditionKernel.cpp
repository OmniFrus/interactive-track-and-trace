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
    const double shoreThreshold = 5000.0; // 5km
    const double minVelocity = 1e-4;
    const double maxVelocity = 3.0;
    const double minCoord = 0.05;

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
            // Closest to West (vertical boundary): normal vel.u = 0, apply xi-formula to vel.v
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

        vel.u = std::clamp(vel.u, -maxVelocity, maxVelocity);
        vel.v = std::clamp(vel.v, -maxVelocity, maxVelocity);

        newLon = longitude + metreToDegrees(vel.u * dt);
        newLat = latitude + metreToDegrees(vel.v * dt);

        newLat = std::clamp(newLat, grid->latMin() + epsilon, grid->latMax() - epsilon);
        newLon = std::clamp(newLon, grid->lonMin() + epsilon, grid->lonMax() - epsilon);
    }

    // Reject motion into coast
    double newShoreDist = grid->getShoreDistance(newLat, newLon);
    if (newShoreDist < 200.0) {  // Within 200m of land
        return {latitude, longitude};  // Stay in current location
    }

    // Optional: fallback if nearly stuck
    double vel2 = vel.u * vel.u + vel.v * vel.v;
    if (vel2 < minVelocity * minVelocity) {
        if (vel2 > 0) {
            double norm = std::sqrt(vel2);
            newLat += metreToDegrees((vel.v / norm) * epsilon);
            newLon += metreToDegrees((vel.u / norm) * epsilon);
        } else {
            newLat += epsilon;
            newLon += epsilon;
        }
    }

    return {newLat, newLon};
}
