
import h5py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load particle trajectory
df = pd.read_csv("single_particle_trajectory.csv")

# Load computation grid
with h5py.File("../../data/grid.h5", "r") as f:
    comp_lat = f["latitude"][:]
    comp_lon = f["longitude"][:]
comp_lon_grid, comp_lat_grid = np.meshgrid(comp_lon, comp_lat)

# Load shoreline grid
with h5py.File("../../data/shore_distance.h5", "r") as f:
    shore_mask = f["mask"][:]
    shore_lat = f["lat"][:]
    shore_lon = f["lon"][:]
shore_mask_2d = shore_mask.reshape(len(shore_lat), len(shore_lon))

# Downsample shoreline grid
stride = 10
shore_lat_ds = shore_lat[::stride]
shore_lon_ds = shore_lon[::stride]
shore_mask_ds = shore_mask_2d[::stride, ::stride]
shore_lon_grid_ds, shore_lat_grid_ds = np.meshgrid(shore_lon_ds, shore_lat_ds)

# Plot everything
plt.figure(figsize=(12, 10))
plt.contourf(shore_lon_grid_ds, shore_lat_grid_ds, shore_mask_ds, levels=[0.5, 1], colors='lightgray', alpha=0.5)
plt.plot(comp_lon_grid, comp_lat_grid, color='lightblue', linewidth=0.5)
plt.plot(comp_lon_grid.T, comp_lat_grid.T, color='lightblue', linewidth=0.5)
plt.plot(df["Longitude"], df["Latitude"], marker='o', markersize=3, label="Particle Trajectory", color='red')
plt.scatter(df["Longitude"].iloc[0], df["Latitude"].iloc[0], color='green', label='Start', zorder=5)
plt.scatter(df["Longitude"].iloc[-1], df["Latitude"].iloc[-1], color='black', label='End', zorder=5)

plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.title("Particle Trajectory vs Computation Grid and Shoreline Mask")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("trajectory_vs_grid_and_shoreline.png", dpi=300)
plt.show()
