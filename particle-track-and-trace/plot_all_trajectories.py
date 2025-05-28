import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# Read the trajectory data
df = pd.read_csv('all_particles_trajectory.csv')

# Create a figure with two subplots
fig = plt.figure(figsize=(20, 10))

# First subplot: Map with all trajectories
ax1 = plt.subplot(121, projection=ccrs.PlateCarree())

# Add map features
ax1.add_feature(cfeature.LAND)
ax1.add_feature(cfeature.OCEAN)
ax1.add_feature(cfeature.COASTLINE)
ax1.add_feature(cfeature.BORDERS, linestyle=':')

# Get unique particle IDs
particle_ids = df['ParticleID'].unique()

# Create a colormap for different particles
colors = plt.cm.rainbow(np.linspace(0, 1, len(particle_ids)))

# Plot each particle's trajectory
for particle_id, color in zip(particle_ids, colors):
    particle_data = df[df['ParticleID'] == particle_id]
    
    # Plot the trajectory
    ax1.plot(particle_data['Longitude'], particle_data['Latitude'], 
             color=color, alpha=0.5, linewidth=1,
             transform=ccrs.PlateCarree())
    
    # Add start and end points
    ax1.scatter(particle_data['Longitude'].iloc[0], particle_data['Latitude'].iloc[0],
                color=color, s=50, transform=ccrs.PlateCarree(), marker='o')
    ax1.scatter(particle_data['Longitude'].iloc[-1], particle_data['Latitude'].iloc[-1],
                color=color, s=50, transform=ccrs.PlateCarree(), marker='x')

# Set the map extent to show all trajectories with some padding
lon_min, lon_max = df['Longitude'].min(), df['Longitude'].max()
lat_min, lat_max = df['Latitude'].min(), df['Latitude'].max()
lon_padding = (lon_max - lon_min) * 0.1
lat_padding = (lat_max - lat_min) * 0.1
ax1.set_extent([lon_min - lon_padding, lon_max + lon_padding,
               lat_min - lat_padding, lat_max + lat_padding])

# Add gridlines
ax1.gridlines(draw_labels=True)

# Add title and legend
ax1.set_title('All Particle Trajectories')

# Second subplot: Distance to shore over time for all particles
ax2 = plt.subplot(122)

# Plot distance to shore for each particle
for particle_id, color in zip(particle_ids, colors):
    particle_data = df[df['ParticleID'] == particle_id]
    ax2.plot(particle_data['Step'], particle_data['DistanceToShore'], 
             color=color, alpha=0.5, linewidth=1,
             label=f'Particle {particle_id}')

ax2.set_xlabel('Time Step')
ax2.set_ylabel('Distance to Shore (meters)')
ax2.set_title('Distance to Shore Over Time for All Particles')
ax2.grid(True)

# Add a legend with a reasonable number of entries
if len(particle_ids) > 20:
    # If there are too many particles, only show some in the legend
    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend(handles[::len(particle_ids)//20], labels[::len(particle_ids)//20],
              title='Particle IDs', loc='upper right', bbox_to_anchor=(1.3, 1))
else:
    ax2.legend(title='Particle IDs', loc='upper right', bbox_to_anchor=(1.3, 1))

# Adjust layout and save
plt.tight_layout()
plt.savefig('all_particles_trajectory_map.png', dpi=300, bbox_inches='tight')
plt.close()

print("Map visualization saved as 'all_particles_trajectory_map.png'")

# Create additional statistics
print("\nTrajectory Statistics:")
print(f"Total number of particles: {len(particle_ids)}")
print(f"Total simulation steps: {df['Step'].max() + 1}")

# Calculate some statistics for each particle
stats = []
for particle_id in particle_ids:
    particle_data = df[df['ParticleID'] == particle_id]
    total_distance = np.sum(np.sqrt(
        np.diff(particle_data['Longitude'])**2 + 
        np.diff(particle_data['Latitude'])**2
    ))
    avg_velocity = np.mean(np.sqrt(
        particle_data['VelocityU']**2 + 
        particle_data['VelocityV']**2
    ))
    min_shore_dist = particle_data['DistanceToShore'].min()
    
    stats.append({
        'ParticleID': particle_id,
        'Total Distance': total_distance,
        'Average Velocity': avg_velocity,
        'Min Distance to Shore': min_shore_dist
    })

# Convert to DataFrame and save statistics
stats_df = pd.DataFrame(stats)
stats_df.to_csv('particle_statistics.csv', index=False)
print("\nDetailed statistics saved to 'particle_statistics.csv'") 