#include <ranges>

#include "UVGrid.h"
#include "readdata.h"
#include "interpolate.h"

using namespace std;

UVGrid::UVGrid(string path) {
  cout << "Reading hydrodynamic_U..." << endl;
  auto us = readHydrodynamicU(path);
  cout << "Reading hydrodynamic_V..." << endl;
  auto vs = readHydrodynamicV(path);
  cout << "Reading grid..." << endl;
  if (us.size() != vs.size()) {
    throw domain_error(sizeError2);
  }

  // Read grid data which now includes shore distances
  auto [times, lats, lons, shoreDistances] = readGrid(path);
  this->times = times;
  this->lats = lats;
  this->lons = lons;
  this->shoreDistances = shoreDistances;

  timeSize = times.size();
  latSize = lats.size();
  lonSize = lons.size();

  size_t gridSize = timeSize * latSize * lonSize;
  if (gridSize != us.size()) {
    throw domain_error(sizeError);
  }

  uvData.reserve(gridSize);

  for (auto vel: views::zip(us, vs)) {
    uvData.push_back(Vel(vel));
  }
}

const Vel &UVGrid::getVelocity(size_t timeIndex, size_t latIndex, size_t lonIndex) const {
  if (timeIndex < 0 || timeIndex >= timeSize
      || latIndex < 0 || latIndex >= latSize
      || lonIndex < 0 || lonIndex >= lonSize) {
    throw std::out_of_range(indexOutOfBounds);
  }
  size_t index = timeIndex * (latSize * lonSize) + latIndex * lonSize + lonIndex;
  return uvData[index];
}

double UVGrid::lonStep() const {
  return lons[1] - lons[0];
}

double UVGrid::latStep() const {
  return lats[1] - lats[0];
}

int UVGrid::timeStep() const {
  return times[1] - times[0];
}

double UVGrid::lonMin() const {
  return this->lons.front();
}

double UVGrid::lonMax() const {
  return this->lons.back();
}

double UVGrid::latMin() const {
  return this->lats.front();
}

double UVGrid::latMax() const {
  return this->lats.back();
}

double UVGrid::timeMax() const {
  return this->times.back();
}

void UVGrid::streamSlice(ostream &os, size_t t) {
  for (int x = 0; x < latSize; x++) {
    for (int y = 0; y < lonSize; y++) {
      auto vel = getVelocity(t, x, y);
      os << vel << " ";
    }
    os << endl;
  }
}

double UVGrid::getShoreDistance(double lat, double lon) const {
  // Find the grid cell containing the point
  int latIndex = (lat - lats[0]) / latStep();
  int lonIndex = (lon - lons[0]) / lonStep();
  
  // Clamp indices to valid range
  latIndex = std::clamp(latIndex, 0, static_cast<int>(latSize - 1));
  lonIndex = std::clamp(lonIndex, 0, static_cast<int>(lonSize - 1));
  
  // Get the shore distance at this grid point
  return shoreDistances[latIndex * lonSize + lonIndex];
}

bool UVGrid::isNearShore(double lat, double lon, double threshold) const {
  return getShoreDistance(lat, lon) <= threshold;
}

