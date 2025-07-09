import h5py

with h5py.File("../../data/grid.h5", "r") as f:
    lats = f["latitude"][:]
    lons = f["longitude"][:]
    times = f["times"][:]

print("Grid Information:")
print(f"Latitude range: {lats[0]} to {lats[-1]}")
print(f"Longitude range: {lons[0]} to {lons[-1]}")
print(f"Time range: {times[0]} to {times[-1]}")
print(f"Latitude step: {lats[1] - lats[0]}")
print(f"Longitude step: {lons[1] - lons[0]}")
print(f"Time step: {times[1] - times[0]}")
print(f"Number of latitude points: {len(lats)}")
print(f"Number of longitude points: {len(lons)}")
print(f"Number of time points: {len(times)}") 