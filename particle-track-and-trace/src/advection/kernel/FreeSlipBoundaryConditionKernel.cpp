#include "FreeSlipBoundaryConditionKernel.h"
#include "../interpolate.h"
#include <cmath>
#include <algorithm>

// Implements a free-slip boundary condition where particles can move freely along the coast
// but cannot move into the land. This is achieved by:
// 1. Projecting the velocity onto the shoreline tangent
// 2. Only allowing movement along the tangent direction
// 3. Using a 2km buffer zone around the coast

FreeSlipBoundaryConditionKernel::FreeSlipBoundaryConditionKernel(
    std::unique_ptr<AdvectionKernel> baseKernel,
    std::shared_ptr<UVGrid> grid)
    : baseKernel(std::move(baseKernel))
    , grid(std::move(grid)) {}

std::pair<double, double> FreeSlipBoundaryConditionKernel::advect(int time, double latitude, double longitude, int dt) const {
    auto vel = bilinearinterpolate(*grid, time, latitude, longitude);

    const double epsilon = 1e-5;
    const double bufferDistance = 5000.0;
    const double delta = AdvectionKernel::metreToDegrees(100); // ~100m

    // Distance to shore — used only to apply BC near coast
    double shoreDist = grid->getShoreDistance(latitude, longitude);

    // Apply BC only inside coastal buffer zone
    if (shoreDist < bufferDistance) {
        // === Compute gradient of speed ===
        auto vel_xp = bilinearinterpolate(*grid, time, latitude, longitude + delta);
        auto vel_xm = bilinearinterpolate(*grid, time, latitude, longitude - delta);
        auto vel_yp = bilinearinterpolate(*grid, time, latitude + delta, longitude);
        auto vel_ym = bilinearinterpolate(*grid, time, latitude - delta, longitude);

        double speed_here = std::sqrt(vel.u * vel.u + vel.v * vel.v);
        double speed_xp = std::sqrt(vel_xp.u * vel_xp.u + vel_xp.v * vel_xp.v);
        double speed_xm = std::sqrt(vel_xm.u * vel_xm.u + vel_xm.v * vel_xm.v);
        double speed_yp = std::sqrt(vel_yp.u * vel_yp.u + vel_yp.v * vel_yp.v);
        double speed_ym = std::sqrt(vel_ym.u * vel_ym.u + vel_ym.v * vel_ym.v);

        double grad_x = (speed_xp - speed_xm) / (2 * delta);
        double grad_y = (speed_yp - speed_ym) / (2 * delta);

        double norm = std::sqrt(grad_x * grad_x + grad_y * grad_y);

        if (norm >= 1e-8) {
            double nX = grad_x / norm;
            double nY = grad_y / norm;
            double tX = -nY;
            double tY = nX;

            // Project velocity onto tangent
            double v_t = vel.u * tX + vel.v * tY;

            double slip_u = v_t * tX;
            double slip_v = v_t * tY;

            double newLon = longitude + AdvectionKernel::metreToDegrees(slip_u * dt);
            double newLat = latitude + AdvectionKernel::metreToDegrees(slip_v * dt);

            newLat = std::clamp(newLat, grid->latMin() + epsilon, grid->latMax() - epsilon);
            newLon = std::clamp(newLon, grid->lonMin() + epsilon, grid->lonMax() - epsilon);

            return {newLat, newLon};
        } else {
            // If gradient too small — prevent movement
            return {latitude, longitude};
        }
    }

    // === Outside buffer zone: normal advection ===
    return baseKernel->advect(time, latitude, longitude, dt);
}