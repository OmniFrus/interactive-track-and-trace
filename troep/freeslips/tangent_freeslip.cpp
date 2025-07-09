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
    // Get the velocity at the current position using bilinear interpolation
    auto vel = bilinearinterpolate(*grid, time, latitude, longitude);

    // Constants for boundary condition calculations
    const double epsilon = 1e-5;          // Small value to prevent particles from exactly hitting boundaries
    const double shoreThreshold = 2000.0; // 2km buffer zone around the coast
    const double minVelocity = 1e-4;      // Minimum velocity threshold
    const double maxVelocity = 3.0;       // Maximum velocity threshold
    const double minCoord = 0.05;         // Minimum coordinate value

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
            
            // Project velocity onto tangent direction only (free-slip)
            double v_t = vel.u * tangent_x + vel.v * tangent_y;
            
            // Calculate slip velocity (only tangential component)
            double slip_u = v_t * tangent_x;
            double slip_v = v_t * tangent_y;
            
            // Update position using slip velocity
            newLon = longitude + metreToDegrees(slip_u * dt);
            newLat = latitude + metreToDegrees(slip_v * dt);
            
            // Ensure particle stays within grid boundaries
            newLat = std::clamp(newLat, grid->latMin() + epsilon, grid->latMax() - epsilon);
            newLon = std::clamp(newLon, grid->lonMin() + epsilon, grid->lonMax() - epsilon);

            return {newLat, newLon};
        } else {
            // If gradient is too small, prevent movement
            return {latitude, longitude};
        }
    }

    // If not near shore, use normal advection
    auto advected = baseKernel->advect(time, latitude, longitude, dt);
    return advected;
}