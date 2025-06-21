import h5py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# === Load shoreline mask ===
with h5py.File("../../data/shore_distance.h5", "r") as f:
    shore_lat = f["lat"][:]
    shore_lon = f["lon"][:]
    shore_mask = f["mask"][:].reshape(len(shore_lat), len(shore_lon))

shore_lon_grid, shore_lat_grid = np.meshgrid(shore_lon, shore_lat)

# === Load velocity grid ===
with h5py.File("../../data/grid.h5", "r") as f:
    grid_lat = f["latitude"][:]
    grid_lon = f["longitude"][:]

grid_lon_grid, grid_lat_grid = np.meshgrid(grid_lon, grid_lat)

# === Load particle trajectory ===
df = pd.read_csv('single_particle_trajectory.csv')  # or all_particles_trajectory.csv

# === Plot ===
fig = plt.figure(figsize=(12, 10))
ax = plt.axes(projection=ccrs.PlateCarree())

# Plot shoreline mask as background
ax.contourf(shore_lon_grid, shore_lat_grid, shore_mask, levels=[0.5, 1], colors='lightgray', alpha=0.8, transform=ccrs.PlateCarree())

# Plot velocity grid lines
ax.plot(grid_lon_grid, grid_lat_grid, color='lightblue', linewidth=0.8, transform=ccrs.PlateCarree())
ax.plot(grid_lon_grid.T, grid_lat_grid.T, color='lightblue', linewidth=0.8, transform=ccrs.PlateCarree())

# Plot particle trajectory
ax.plot(df['Longitude'], df['Latitude'], 'r-', linewidth=2, label='Trajectory', transform=ccrs.PlateCarree())

# Start and end points
ax.scatter(df['Longitude'].iloc[0], df['Latitude'].iloc[0], color='green', s=100, label='Start', transform=ccrs.PlateCarree())
ax.scatter(df['Longitude'].iloc[-1], df['Latitude'].iloc[-1], color='black', s=100, label='End', transform=ccrs.PlateCarree())

# Set extent to grid boundaries
# Compute particle bounding box with padding
lon_min, lon_max = df['Longitude'].min(), df['Longitude'].max()
lat_min, lat_max = df['Latitude'].min(), df['Latitude'].max()

lon_padding = (lon_max - lon_min) * 0.1
lat_padding = (lat_max - lat_min) * 0.1

ax.set_extent([lon_min - lon_padding, lon_max + lon_padding,
               lat_min - lat_padding, lat_max + lat_padding])

# Title and legend
ax.set_title('Particle Trajectory with Shoreline Mask and Velocity Grid')
ax.legend()

# Save figure
plt.savefig('trajectory_with_shoreline_and_grid.png', dpi=300, bbox_inches='tight')
plt.show()
