#include <stdexcept>

#include <netcdf>

#include "readdata.h"

using namespace std;
using namespace netCDF;

template<typename T>
vector<T> getVarVector(const NcVar &var) {
  int length = 1;
  for (NcDim dim: var.getDims()) {
    length *= dim.getSize();
  }

  vector<T> vec(length);

  var.getVar(vec.data());

  return vec;
}

vector<double> readHydrodynamicU(string path) {
  // Vs and Us flipped cause the files are named incorrectly
  string fileName = "hydrodynamic_V.h5";
  cout << "Opening " << path + '/' + fileName << endl;
  netCDF::NcFile data(path + '/' + fileName, netCDF::NcFile::read);

  multimap<string, NcVar> vars = data.getVars();

  return getVarVector<double>(vars.find("vo")->second);
}

vector<double> readHydrodynamicV(string path) {
  // Vs and Us flipped cause the files are named incorrectly
  string fileName = "hydrodynamic_U.h5";
  cout << "Opening " << path + '/' + fileName << endl;
  netCDF::NcFile data(path + '/' + fileName, netCDF::NcFile::read);

  multimap<string, NcVar> vars = data.getVars();

  return getVarVector<double>(vars.find("uo")->second);
}

tuple<vector<int>, vector<double>, vector<double>> readGrid(string path) {
  string fileName = "grid.h5";
  netCDF::NcFile data(path + '/' + fileName, netCDF::NcFile::read);
  multimap<string, NcVar> vars = data.getVars();
  vector<int> time = getVarVector<int>(vars.find("times")->second);
  vector<double> longitude = getVarVector<double>(vars.find("longitude")->second);
  vector<double> latitude = getVarVector<double>(vars.find("latitude")->second);

  return {time, latitude, longitude};
}

tuple<vector<double>, vector<double>, vector<double>, vector<uint8_t>> readShoreDistance(string path) {
  string fileName = "shore_distance.h5";
  cout << "Opening " << path + '/' + fileName << endl;
  netCDF::NcFile data(path + '/' + fileName, netCDF::NcFile::read);
  multimap<string, NcVar> vars = data.getVars();

  auto it = vars.find("distance");
  if (it == vars.end()) {
      throw std::runtime_error("'distance' not found in shore_distance.h5");
  }
  vector<double> distances = getVarVector<double>(it->second);
  vector<double> latitude = getVarVector<double>(vars.find("lat")->second);
  vector<double> longitude = getVarVector<double>(vars.find("lon")->second);
  vector<uint8_t> mask = getVarVector<uint8_t>(vars.find("mask")->second);

  return {distances, latitude, longitude, mask};
}
