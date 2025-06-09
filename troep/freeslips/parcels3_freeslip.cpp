#include "parcels3_freeslip.h"
#include "../interpolate.h"
#include <cmath>
#include <algorithm>

ParcelsRK4FreeSlipKernel::ParcelsRK4FreeSlipKernel(std::shared_ptr<UVGrid> grid)
    : grid(std::move(grid)) {}

std::pair<double, double> ParcelsRK4FreeSlipKernel::advect(int time, double lat, double lon, int dt) const {
    const double epsilon = 1e-5;
    const double maxVelocity = 3.0;

    auto applySlip = [&](double u, double v, double lat, double lon) -> Vel {
        double dLat = grid->latStep();
        double dLon = grid->lonStep();

        double xi  = (lon - grid->lonMin()) / dLon;
        double eta = (lat - grid->latMin()) / dLat;

        xi  = xi - std::floor(xi);
        eta = eta - std::floor(eta);

        xi  = std::clamp(xi, epsilon, 1.0 - epsilon);
        eta = std::clamp(eta, epsilon, 1.0 - epsilon);

        bool landE = grid->isLand(lat, lon + dLon);
        bool landW = grid->isLand(lat, lon - dLon);
        bool landN = grid->isLand(lat + dLat, lon);
        bool landS = grid->isLand(lat - dLat, lon);

        if (landW && !landE) { u = 0.0; v *= 1.0 / xi; }
        if (landE && !landW) { u = 0.0; v *= 1.0 / (1.0 - xi); }
        if (landS && !landN) { v = 0.0; u *= 1.0 / eta; }
        if (landN && !landS) { v = 0.0; u *= 1.0 / (1.0 - eta); }

        u = std::clamp(u, -maxVelocity, maxVelocity);
        v = std::clamp(v, -maxVelocity, maxVelocity);

        return {u, v};
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
