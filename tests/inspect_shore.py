import h5py
import numpy as np
import matplotlib.pyplot as plt

# Path to your shore_distance.h5 file
filename = "data/shore_distance.h5"

with h5py.File(filename, "r") as f:
    shore_dist = f["distance"][:]
    shore_lat = f["lat"][:]
    shore_lon = f["lon"][:]
    land_mask = f["mask"][:].reshape(shore_dist.shape)

lon_grid, lat_grid = np.meshgrid(shore_lon, shore_lat)

plt.figure(figsize=(10, 6))
pcm = plt.pcolormesh(lon_grid, lat_grid, shore_dist, shading="auto", cmap="viridis")
plt.colorbar(pcm, label="Distance to Shore (m)")
plt.contour(lon_grid, lat_grid, land_mask, levels=[0.5], colors='red', linewidths=1.5)
plt.title("Shore Distance with Land Mask Overlay (red = land)")
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.tight_layout()
plt.show()
