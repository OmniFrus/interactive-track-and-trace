import h5py

with h5py.File("data/shore_distance.h5", "r") as f:
    mask = f["mask"][:]
    lats = f["lat"][:]
    lons = f["lon"][:]

print("Land Mask Information:")
print(f"Land mask shape: {mask.shape}")
print(f"Number of land cells: {mask.sum()}")
print(f"Number of water cells: {mask.size - mask.sum()}")
print(f"Latitude range: {lats[0]} to {lats[-1]}")
print(f"Longitude range: {lons[0]} to {lons[-1]}")
print(f"Latitude step: {lats[1] - lats[0]}")
print(f"Longitude step: {lons[1] - lons[0]}")
print(f"Number of latitude points: {len(lats)}")
print(f"Number of longitude points: {len(lons)}") 