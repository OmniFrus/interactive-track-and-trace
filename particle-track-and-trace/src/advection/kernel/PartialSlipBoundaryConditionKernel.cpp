#include "PartialSlipBoundaryConditionKernel.h"
#include "../interpolate.h"
#include <cmath>
#include <algorithm>

// Implements a partial-slip boundary condition where particles can move both along the coast
// and partially into/out of the land. This is achieved by:
// 1. Projecting the velocity onto normal and tangential components
// 2. Reducing the normal component by a slip ratio (0.5)
// 3. Keeping the full tangential component
// 4. Using a 2km buffer zone around the coast

PartialSlipBoundaryConditionKernel::PartialSlipBoundaryConditionKernel(
    std::unique_ptr<AdvectionKernel> baseKernel,
    std::shared_ptr<UVGrid> grid)
    : baseKernel(std::move(baseKernel))
    , grid(std::move(grid)) {}

std::pair<double, double> PartialSlipBoundaryConditionKernel::advect(int time, double latitude, double longitude, int dt) const {
    // Get the velocity at the current position using bilinear interpolation
    auto vel = bilinearinterpolate(*grid, time, latitude, longitude);

    // Constants for boundary condition calculations
    const double epsilon = 1e-5;          // Small value to prevent particles from exactly hitting boundaries
    const double shoreThreshold = 2000.0; // 2km buffer zone around the coast
    const double minVelocity = 1e-4;      // Minimum velocity threshold
    const double maxVelocity = 3.0;       // Maximum velocity threshold
    const double slipRatio = 0.5;         // Moderate slip (0 = no slip, 1 = full slip)

    // Initialize new position with current position
    double newLat = latitude;
    double newLon = longitude;
    
    // Get distance to shore at current position
    double shoreDist = grid->getShoreDistance(latitude, longitude);
    
    // Only apply boundary conditions if we're within the shore buffer zone
    if (shoreDist < shoreThreshold) {
        // Calculate gradient of shore distance to determine normal and tangential directions
        double delta = 1e-4; // Small increment in degrees for numerical gradient
        double gradLat = (grid->getShoreDistance(latitude + delta, longitude) - 
                         grid->getShoreDistance(latitude - delta, longitude)) / (2 * delta);
        double gradLon = (grid->getShoreDistance(latitude, longitude + delta) - 
                         grid->getShoreDistance(latitude, longitude - delta)) / (2 * delta);
        
        // Calculate norm of gradient vector for normalization
        double norm = std::sqrt(gradLat * gradLat + gradLon * gradLon);
        
        // Only proceed if gradient is significant enough
        if (norm >= 1e-8) {
            // Calculate unit normal vector (pointing away from shore)
            double normal_x = gradLon / norm;
            double normal_y = gradLat / norm;
            
            // Calculate unit tangent vector (perpendicular to normal)
            double tangent_x = -normal_y;
            double tangent_y = normal_x;
            
            // Project velocity onto normal and tangential components
            double v_n = vel.u * normal_x + vel.v * normal_y;  // Normal component
            double v_t = vel.u * tangent_x + vel.v * tangent_y; // Tangential component

            // Calculate slip velocity:
            // Normal component is reduced by (1-slipRatio)
            // Tangential component remains full
            double slip_u = (1.0 - slipRatio) * v_n * normal_x + v_t * tangent_x;
            double slip_v = (1.0 - slipRatio) * v_n * normal_y + v_t * tangent_y;
            
            // Update position using slip velocity
            newLon = longitude + metreToDegrees(slip_u * dt);
            newLat = latitude + metreToDegrees(slip_v * dt);
            
            // Ensure particle stays within grid boundaries
            newLat = std::clamp(newLat, grid->latMin() + epsilon, grid->latMax() - epsilon);
            newLon = std::clamp(newLon, grid->lonMin() + epsilon, grid->lonMax() - epsilon);

            return {newLat, newLon};
        }
    }

    // If not near shore, use normal advection
    auto advected = baseKernel->advect(time, latitude, longitude, dt);
    return advected;
}