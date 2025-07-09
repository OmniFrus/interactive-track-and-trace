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
    shore_distance = f["distance"][:].reshape(len(shore_lat), len(shore_lon))

shore_lon_grid, shore_lat_grid = np.meshgrid(shore_lon, shore_lat)

# === Load velocity grid ===
with h5py.File("../../data/grid.h5", "r") as f:
    grid_lat = f["latitude"][:]
    grid_lon = f["longitude"][:]

grid_lon_grid, grid_lat_grid = np.meshgrid(grid_lon, grid_lat)

# === Load particle trajectory ===
df = pd.read_csv('single_particle_trajectory.csv')

# === Create 2 subplots ===
fig = plt.figure(figsize=(15, 10))

# First subplot: Map
ax1 = plt.subplot(121, projection=ccrs.PlateCarree())

# Zoom on particle
lon_min, lon_max = df['Longitude'].min(), df['Longitude'].max()
lat_min, lat_max = df['Latitude'].min(), df['Latitude'].max()
lon_padding = (lon_max - lon_min) * 0.1
lat_padding = (lat_max - lat_min) * 0.1

ax1.set_extent([lon_min - lon_padding, lon_max + lon_padding,
                lat_min - lat_padding, lat_max + lat_padding])

# Crop shore_distance grid using same padding:
bbox = [lon_min - lon_padding, lon_max + lon_padding,
        lat_min - lat_padding, lat_max + lat_padding]

lon_mask = (shore_lon >= bbox[0]) & (shore_lon <= bbox[1])
lat_mask = (shore_lat >= bbox[2]) & (shore_lat <= bbox[3])

shore_lon_crop = shore_lon[lon_mask]
shore_lat_crop = shore_lat[lat_mask]
shore_distance_crop = shore_distance[np.ix_(lat_mask, lon_mask)]

shore_lon_grid_crop, shore_lat_grid_crop = np.meshgrid(shore_lon_crop, shore_lat_crop)

# Now plot cropped GEBCO grid:
c = ax1.pcolormesh(shore_lon_grid_crop, shore_lat_grid_crop, shore_distance_crop,
                   cmap='Blues_r', shading='auto', transform=ccrs.PlateCarree())

plt.colorbar(c, ax=ax1, label='Distance to Shore (meters)')

# Plot shoreline mask
ax1.contourf(shore_lon_grid, shore_lat_grid, shore_mask, levels=[0.5, 1], colors='lightgray', alpha=0.8, transform=ccrs.PlateCarree())

# Velocity grid lines
ax1.plot(grid_lon_grid, grid_lat_grid, color='lightblue', linewidth=0.8, transform=ccrs.PlateCarree())
ax1.plot(grid_lon_grid.T, grid_lat_grid.T, color='lightblue', linewidth=0.8, transform=ccrs.PlateCarree())

# Particle trajectory
ax1.plot(df['Longitude'], df['Latitude'], 'r-', linewidth=2, label='Trajectory', transform=ccrs.PlateCarree())

# Start and end points
ax1.scatter(df['Longitude'].iloc[0], df['Latitude'].iloc[0], color='green', s=100, label='Start', transform=ccrs.PlateCarree())
ax1.scatter(df['Longitude'].iloc[-1], df['Latitude'].iloc[-1], color='black', s=100, label='End', transform=ccrs.PlateCarree())


ax1.set_title('Particle Trajectory with Shoreline Mask and Velocity Grid')
ax1.legend()

# Second subplot: Distance to Shore
ax2 = plt.subplot(122)
ax2.plot(df['Step'], df['DistanceToShore'], 'b-', label='Distance to Shore')
ax2.set_xlabel('Time Step')
ax2.set_ylabel('Distance to Shore (meters)')
ax2.set_title('Distance to Shore Over Time')
ax2.grid(True)

if df['Beached'].eq('Yes').any():
    beach_step = df[df['Beached'] == 'Yes']['Step'].iloc[0]
    ax2.axvline(x=beach_step, linestyle='--', color='purple', label=f'Beached at Step {beach_step}')

ax2.legend()

# Finalize
plt.tight_layout()
plt.savefig('particle_trajectory_map.png', dpi=300, bbox_inches='tight')
plt.show()