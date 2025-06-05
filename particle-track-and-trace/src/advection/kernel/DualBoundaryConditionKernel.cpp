#include "DualBoundaryConditionKernel.h"
#include "../interpolate.h"
#include <cmath>
#include <algorithm>

// Basically a partial slip with a slip factor close to 1, but with different beaching mechanics.
// These beaching mechanics are in LagrangeGlyphs.cpp

DualBoundaryConditionKernel::DualBoundaryConditionKernel(
    std::unique_ptr<AdvectionKernel> baseKernel,
    std::shared_ptr<UVGrid> grid)
    : baseKernel(std::move(baseKernel))
    , grid(std::move(grid)) {}

std::pair<double, double> DualBoundaryConditionKernel::advect(int time, double latitude, double longitude, int dt) const {
    auto vel = bilinearinterpolate(*grid, time, latitude, longitude);

    const double epsilon = 1e-5;
    const double shoreThreshold = 2000.0; // 2km
    const double minVelocity = 1e-4;
    const double maxVelocity = 3.0;
    const double slipRatio = 0.95;  // Between 0 (no slip) and 1 (full slip)

    // Only apply partial slip: blend between normal and tangential movement if near shore
    double newLat = latitude;
    double newLon = longitude;
    double shoreDist = grid->getShoreDistance(latitude, longitude);
    if (shoreDist < shoreThreshold) { // Only apply if near shore (buffer)
        double delta = 1e-4; // Small increment in degrees
        double gradLat = (grid->getShoreDistance(latitude + delta, longitude) - grid->getShoreDistance(latitude - delta, longitude)) / (2 * delta);
        double gradLon = (grid->getShoreDistance(latitude, longitude + delta) - grid->getShoreDistance(latitude, longitude - delta)) / (2 * delta);
        double norm = std::sqrt(gradLat * gradLat + gradLon * gradLon);
        if (norm >= 1e-8) {
            double normal_x = gradLon / norm;
            double normal_y = gradLat / norm;
            double tangent_x = -normal_y;
            double tangent_y = normal_x;
            
            double v_n = vel.u * normal_x + vel.v * normal_y;
            double v_t = vel.u * tangent_x + vel.v * tangent_y;

            // Keep full tangential component, reduce normal component based on slip ratio
            double slip_u = (1.0 - slipRatio) * v_n * normal_x + v_t * tangent_x;
            double slip_v = (1.0 - slipRatio) * v_n * normal_y + v_t * tangent_y;
            
            newLon = longitude + metreToDegrees(slip_u * dt);
            newLat = latitude + metreToDegrees(slip_v * dt);
            newLat = std::clamp(newLat, grid->latMin() + epsilon, grid->latMax() - epsilon);
            newLon = std::clamp(newLon, grid->lonMin() + epsilon, grid->lonMax() - epsilon);

            return {newLat, newLon};
        }
    }

    // If not near shore, normal advection
    auto advected = baseKernel->advect(time, latitude, longitude, dt);
    return advected;
}
