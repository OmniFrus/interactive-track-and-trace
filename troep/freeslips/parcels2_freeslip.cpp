#include "FreeSlipBoundaryConditionKernel.h"
#include "../interpolate.h"
#include <cmath>
#include <algorithm>

FreeSlipBoundaryConditionKernel::FreeSlipBoundaryConditionKernel(
    std::unique_ptr<AdvectionKernel> baseKernel,
    std::shared_ptr<UVGrid> grid)
    : baseKernel(std::move(baseKernel))
    , grid(std::move(grid)) {}

std::pair<double, double> FreeSlipBoundaryConditionKernel::advect(int time, double lat, double lon, int dt) const {
    const double epsilon = 1e-5;
    const double maxVelocity = 3.0;

    // Grid resolution
    double dLat = grid->latStep();
    double dLon = grid->lonStep();

    // Velocity at this location
    auto vel = bilinearinterpolate(*grid, time, lat, lon);

    // Compute fractional index coordinates
    double xi  = (lon - grid->lonMin()) / dLon;
    double eta = (lat - grid->latMin()) / dLat;

    xi  = xi - std::floor(xi);   // fractional part [0,1)
    eta = eta - std::floor(eta); // fractional part [0,1)

    // Avoid division by zero
    xi  = std::clamp(xi, epsilon, 1.0 - epsilon);
    eta = std::clamp(eta, epsilon, 1.0 - epsilon);

    // Land detection around particle
    bool landE = grid->isLand(lat, lon + dLon);
    bool landW = grid->isLand(lat, lon - dLon);
    bool landN = grid->isLand(lat + dLat, lon);
    bool landS = grid->isLand(lat - dLat, lon);

    // Apply Parcels-style free-slip correction (independent in each direction)
    if (landW && !landE) {
        vel.u = 0.0;
        vel.v *= 1.0 / xi;
    }
    if (landE && !landW) {
        vel.u = 0.0;
        vel.v *= 1.0 / (1.0 - xi);
    }
    if (landS && !landN) {
        vel.v = 0.0;
        vel.u *= 1.0 / eta;
    }
    if (landN && !landS) {
        vel.v = 0.0;
        vel.u *= 1.0 / (1.0 - eta);
    }

    // Clamp extreme velocities
    vel.u = std::clamp(vel.u, -maxVelocity, maxVelocity);
    vel.v = std::clamp(vel.v, -maxVelocity, maxVelocity);

    // Integrate position
    double newLon = lon + metreToDegrees(vel.u * dt);
    double newLat = lat + metreToDegrees(vel.v * dt);

    // Clamp to grid bounds
    newLat = std::clamp(newLat, grid->latMin() + epsilon, grid->latMax() - epsilon);
    newLon = std::clamp(newLon, grid->lonMin() + epsilon, grid->lonMax() - epsilon);

    if (grid->isLand(newLat, newLon)) {
        // Compute gradient of shore distance field
        const double delta = metreToDegrees(100.0);  // 100m in degrees

        double gradLat = (grid->getShoreDistance(newLat + delta, newLon) -
                        grid->getShoreDistance(newLat - delta, newLon)) / (2.0 * delta);

        double gradLon = (grid->getShoreDistance(newLat, newLon + delta) -
                        grid->getShoreDistance(newLat, newLon - delta)) / (2.0 * delta);

        // Normalize the gradient
        double norm = std::sqrt(gradLat * gradLat + gradLon * gradLon);
        if (norm > 0.0) {
            gradLat /= norm;
            gradLon /= norm;
        } else {
            gradLat = gradLon = 0.0;
        }

        // Push particle slightly away from shore
        double pushedLat = lat + epsilon * gradLat;
        double pushedLon = lon + epsilon * gradLon;

        // If still in land, just stay in place
        if (!grid->isLand(pushedLat, pushedLon)) {
            return {pushedLat, pushedLon};
        } else {
            return {lat, lon};
        }
    }  
    
    return {newLat, newLon};
}