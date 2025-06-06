import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import h5py

# Read the trajectory data
df = pd.read_csv('single_particle_trajectory.csv')

# Create a figure with two subplots
fig = plt.figure(figsize=(15, 10))

# First subplot: Map with trajectory
ax1 = plt.subplot(121, projection=ccrs.PlateCarree())

# Add map features
ax1.add_feature(cfeature.LAND)
ax1.add_feature(cfeature.OCEAN)
ax1.add_feature(cfeature.COASTLINE)
ax1.add_feature(cfeature.BORDERS, linestyle=':')

# Plot the trajectory with color based on distance to shore (meters)
scatter = ax1.scatter(df['Longitude'], df['Latitude'], 
                     c=df['DistanceToShore'], cmap='viridis',
                     transform=ccrs.PlateCarree(),
                     label='Particle Path')

# Add colorbar
plt.colorbar(scatter, ax=ax1, label='Distance to Shore (meters)')

# Add start and end points
ax1.scatter(df['Longitude'].iloc[0], df['Latitude'].iloc[0], 
           color='green', s=100, transform=ccrs.PlateCarree(),
           label='Start')
ax1.scatter(df['Longitude'].iloc[-1], df['Latitude'].iloc[-1], 
           color='red', s=100, transform=ccrs.PlateCarree(),
           label='End')

# Add velocity vectors (every 10th point to avoid overcrowding)
step = 10
for i in range(0, len(df), step):
    ax1.quiver(df['Longitude'].iloc[i], df['Latitude'].iloc[i],
              df['VelocityU'].iloc[i], df['VelocityV'].iloc[i],
              color='black', scale=50, transform=ccrs.PlateCarree())
    
# Overlay grid points from grid.h5 or shore_distance.h5
grid_file = "../../data/grid.h5"  # or "shore_distance.h5" if you prefer
with h5py.File(grid_file, "r") as f:
    # Try both naming conventions
    if "latitude" in f and "longitude" in f:
        lats = f["latitude"][:]
        lons = f["longitude"][:]
    elif "lat" in f and "lon" in f:
        lats = f["lat"][:]
        lons = f["lon"][:]
    else:
        raise ValueError("No lat/lon found in file")

# Create meshgrid for all grid points
lon_grid, lat_grid = np.meshgrid(lons, lats)

# Overlay grid points on the map
ax1.scatter(lon_grid, lat_grid, s=8, color='black', alpha=0.4, label='Grid Points', transform=ccrs.PlateCarree())


# Overlay grid points from grid.h5 or shore_distance.h5
grid_file = "../../data/grid.h5"  # or "shore_distance.h5" if you prefer
with h5py.File(grid_file, "r") as f:
    # Try both naming conventions
    if "latitude" in f and "longitude" in f:
        lats = f["latitude"][:]
        lons = f["longitude"][:]
    elif "lat" in f and "lon" in f:
        lats = f["lat"][:]
        lons = f["lon"][:]
    else:
        raise ValueError("No lat/lon found in file")

# Create meshgrid for all grid points
lon_grid, lat_grid = np.meshgrid(lons, lats)

# Overlay grid points on the map
ax1.scatter(lon_grid, lat_grid, s=8, color='black', alpha=0.4, label='Grid Points', transform=ccrs.PlateCarree())

# Set the map extent to show the full trajectory with some padding
lon_min, lon_max = df['Longitude'].min(), df['Longitude'].max()
lat_min, lat_max = df['Latitude'].min(), df['Latitude'].max()
lon_padding = (lon_max - lon_min) * 0.1
lat_padding = (lat_max - lat_min) * 0.1
ax1.set_extent([lon_min - lon_padding, lon_max + lon_padding,
               lat_min - lat_padding, lat_max + lat_padding])

# Add gridlines
ax1.gridlines(draw_labels=True)

# Add title and legend
ax1.set_title('Particle Trajectory with Distance to Shore')
ax1.legend()

# Second subplot: Distance to shore over time
ax2 = plt.subplot(122)
ax2.plot(df['Step'], df['DistanceToShore'], 'b-', label='Distance to Shore')
ax2.set_xlabel('Time Step')
ax2.set_ylabel('Distance to Shore (meters)')
ax2.set_title('Distance to Shore Over Time')
ax2.grid(True)
ax2.legend()

# Adjust layout and save
plt.tight_layout()
plt.savefig('particle_trajectory_map.png', dpi=300, bbox_inches='tight')
plt.close()

print("Map visualization saved as 'particle_trajectory_map.png'") 