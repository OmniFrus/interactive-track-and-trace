import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
import h5py

# Read the trajectory data
df = pd.read_csv('all_particles_trajectory.csv')

# === Load shoreline mask ===
with h5py.File("../../data/shore_distance.h5", "r") as f:
    shore_lat = f["lat"][:]
    shore_lon = f["lon"][:]
    shore_mask = f["mask"][:].reshape(len(shore_lat), len(shore_lon))

shore_lon_grid, shore_lat_grid = np.meshgrid(shore_lon, shore_lat)

# Create figure
fig = plt.figure(figsize=(12, 8))

# Create subplot for map
ax1 = plt.subplot(111, projection=ccrs.PlateCarree())

# Plot shoreline mask
ax1.contourf(shore_lon_grid, shore_lat_grid, shore_mask,
             levels=[0.5, 1], colors='lightgray', alpha=0.8, transform=ccrs.PlateCarree())

# Get unique particle IDs
particle_ids = df['ParticleID'].unique()

# Create colormap
colors = plt.cm.rainbow(np.linspace(0, 1, len(particle_ids)))

# Plot trajectories
for particle_id, color in zip(particle_ids, colors):
    particle_data = df[df['ParticleID'] == particle_id]
    
    # Plot trajectory
    ax1.plot(particle_data['Longitude'], particle_data['Latitude'], 
             color=color, alpha=0.5, linewidth=1,
             transform=ccrs.PlateCarree())
    
    # Add start and end markers
    ax1.scatter(particle_data['Longitude'].iloc[0], particle_data['Latitude'].iloc[0],
                color=color, s=50, transform=ccrs.PlateCarree(), marker='o')
    ax1.scatter(particle_data['Longitude'].iloc[-1], particle_data['Latitude'].iloc[-1],
                color=color, s=50, transform=ccrs.PlateCarree(), marker='x')
    
    # Mark beaching point if it exists
    if particle_data['Beached'].eq('Yes').any():
        beach_index = particle_data[particle_data['Beached'] == 'Yes'].index[0]
        beached_step = particle_data.loc[beach_index, 'Step']
        beached_lon = particle_data.loc[beach_index, 'Longitude']
        beached_lat = particle_data.loc[beach_index, 'Latitude']
        ax1.scatter(beached_lon, beached_lat, color=color, s=100, 
                    transform=ccrs.PlateCarree(), marker='*')
        ax1.annotate(f'Step {beached_step}', 
                     (beached_lon, beached_lat),
                     xytext=(10, 10), textcoords='offset points',
                     color=color, fontsize=8)

# Set map extent
lon_min, lon_max = df['Longitude'].min(), df['Longitude'].max()
lat_min, lat_max = df['Latitude'].min(), df['Latitude'].max()
lon_padding = (lon_max - lon_min) * 0.1
lat_padding = (lat_max - lat_min) * 0.1
ax1.set_extent([lon_min - lon_padding, lon_max + lon_padding,
               lat_min - lat_padding, lat_max + lat_padding])

# Add gridlines
ax1.gridlines(draw_labels=True)

# Add title
ax1.set_title('All Particle Trajectories\n* indicates beaching location')

# Save the plot
plt.tight_layout()
plt.savefig('all_particles_trajectory_map.png', dpi=300, bbox_inches='tight')
plt.close()

print("Trajectory map saved as 'all_particles_trajectory_map.png'")
