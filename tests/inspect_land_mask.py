import h5py

with h5py.File("data/shore_distance.h5", "r") as f:
    mask = f["mask"][:]
    lats = f["lat"][:]
    lons = f["lon"][:]

# Particle's position
lat = 57.1603
lon = 6.00393

# Find the closest grid point
lat_index = min(range(len(lats)), key=lambda i: abs(lats[i] - lat))
lon_index = min(range(len(lons)), key=lambda i: abs(lons[i] - lon))

print(f"Particle's position: ({lat}, {lon})")
print(f"Closest grid point: ({lats[lat_index]}, {lons[lon_index]})")
print(f"Land mask value: {mask[lat_index, lon_index]}")

# Print the land mask values around the particle's position
print("\nLand mask values around the particle's position:")
for i in range(lat_index - 2, lat_index + 3):
    for j in range(lon_index - 2, lon_index + 3):
        if 0 <= i < len(lats) and 0 <= j < len(lons):
            print(f"({lats[i]}, {lons[j]}): {mask[i, j]}") 