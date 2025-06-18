#include "ParcelsBoundaryConditionKernel.h"
#include "../interpolate.h"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iomanip>

ParcelsRK4FreeSlipKernel::ParcelsRK4FreeSlipKernel(std::unique_ptr<RK4AdvectionKernel> baseKernel, std::shared_ptr<UVGrid> grid)
    : baseKernel(std::move(baseKernel)), grid(grid) {}

std::pair<double, double> ParcelsRK4FreeSlipKernel::advect(int time, double lat, double lon, int dt) const {
    const double epsilon = 1e-5;
    const double maxVelocity = 3.0;

//    static std::ofstream slip_log;
//    static bool initialized = false;
//    if (!initialized) {
//        slip_log.open("slip_violations.csv", std::ios::app);
//        if (!slip_log.is_open()) { 
//           std::cerr << "❌ Failed to open slip_violations.csv for writing!" << std::endl;
//        } else {
//            slip_log << "Latitude,Longitude,DotProduct,VelocityU,VelocityV,NormalX,NormalY\n";
//        }
//        initialized = true;
//    }

    auto applySlip = [&](double u, double v, double lat, double lon) -> Vel {
        const double delta = AdvectionKernel::metreToDegrees(100); // 100m
        const double slipRatio = 0.5; // Full slip; use <1.0 for partial slip
        const double slipRadius = 5000.0; // Only apply slip if within 2km of shore

        //if (grid->getShoreDistance(lat, lon) > slipRadius) { // -> comment out to have unconditional
        //    return {u, v}; // use original velocity
        //}

        double gradLat = (grid->getShoreDistance(lat + delta, lon) -
                          grid->getShoreDistance(lat - delta, lon)) / (2 * delta);
        double gradLon = (grid->getShoreDistance(lat, lon + delta) -
                          grid->getShoreDistance(lat, lon - delta)) / (2 * delta);

        double norm = std::sqrt(gradLat * gradLat + gradLon * gradLon);
        if (norm < 1e-8) return {u, v};  // skip if flat

        // Normal and tangent
        double nX = gradLon / norm;
        double nY = gradLat / norm;
        double tX = -nY;
        double tY = nX;

        double v_n = u * nX + v * nY;
        double v_t = u * tX + v * tY;

        double slip_u;
        double slip_v;

        //if(v_n >= 0.0) { -> Tried something, particles would circle in place. TODO - find a fix
        //    // moving away from shore - no slip
        //    return{u,v};
        //} else { // moving towards shore - slip
        //    slip_u = v_t * tX + (1.0 - slipRatio) * v_n * nX;
        //    slip_v = v_t * tY + (1.0 - slipRatio) * v_n * nY;
        //}

        slip_u = v_t * tX + (1.0 - slipRatio) * v_n * nX; // comment these when using the if-statement above.
        slip_v = v_t * tY + (1.0 - slipRatio) * v_n * nY;

        // Log violation if slip is not tangential
//        double finalDotWithNormal = slip_u * nX + slip_v * nY;
//        if (std::abs(finalDotWithNormal) > 1e-5) {
//            slip_log << std::fixed << std::setprecision(8)
//                     << lat << "," << lon << ","
//                     << finalDotWithNormal << ","
//                     << u << "," << v << ","
//                     << nX << "," << nY << "\n";
//        }

        return {slip_u, slip_v};
    };

    auto [u1, v1] = applySlip(bilinearinterpolate(*grid, time, lat, lon).u,
                              bilinearinterpolate(*grid, time, lat, lon).v,
                              lat, lon);
    double lon1 = lon + metreToDegrees(u1 * 0.5 * dt);
    double lat1 = lat + metreToDegrees(v1 * 0.5 * dt);

    auto [u2, v2] = applySlip(bilinearinterpolate(*grid, time + 0.5 * dt, lat1, lon1).u,
                              bilinearinterpolate(*grid, time + 0.5 * dt, lat1, lon1).v,
                              lat1, lon1);
    double lon2 = lon + metreToDegrees(u2 * 0.5 * dt);
    double lat2 = lat + metreToDegrees(v2 * 0.5 * dt);

    auto [u3, v3] = applySlip(bilinearinterpolate(*grid, time + 0.5 * dt, lat2, lon2).u,
                              bilinearinterpolate(*grid, time + 0.5 * dt, lat2, lon2).v,
                              lat2, lon2);
    double lon3 = lon + metreToDegrees(u3 * dt);
    double lat3 = lat + metreToDegrees(v3 * dt);

    auto [u4, v4] = applySlip(bilinearinterpolate(*grid, time + dt, lat3, lon3).u,
                              bilinearinterpolate(*grid, time + dt, lat3, lon3).v,
                              lat3, lon3);

    double lonFinal = lon + metreToDegrees((u1 + 2 * u2 + 2 * u3 + u4) / 6.0 * dt);
    double latFinal = lat + metreToDegrees((v1 + 2 * v2 + 2 * v3 + v4) / 6.0 * dt);

    latFinal = std::clamp(latFinal, grid->latMin() + epsilon, grid->latMax() - epsilon);
    lonFinal = std::clamp(lonFinal, grid->lonMin() + epsilon, grid->lonMax() - epsilon);

    return {latFinal, lonFinal};
}