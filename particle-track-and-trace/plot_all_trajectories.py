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
    
    # Check if particle is beached (distance to shore is very small)
    beached_mask = particle_data['DistanceToShore'] < 0.001  # 1 meter threshold
    if beached_mask.any():
        beached_step = particle_data.loc[beached_mask.idxmax(), 'Step']
        beached_lon = particle_data.loc[beached_mask.idxmax(), 'Longitude']
        beached_lat = particle_data.loc[beached_mask.idxmax(), 'Latitude']
        ax1.scatter(beached_lon, beached_lat, color=color, s=100, 
                   transform=ccrs.PlateCarree(), marker='*')
        ax1.annotate(f'Step {beached_step}', 
                    (beached_lon, beached_lat),
                    xytext=(10, 10), textcoords='offset points',
                    color=color, fontsize=8)

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
ax1.set_title('All Particle Trajectories\n* indicates beaching location')

# Second subplot: Distance to shore over time for all particles
ax2 = plt.subplot(122)

# Plot distance to shore for each particle
for particle_id, color in zip(particle_ids, colors):
    particle_data = df[df['ParticleID'] == particle_id]
    ax2.plot(particle_data['Step'], particle_data['DistanceToShore'], 
             color=color, alpha=0.5, linewidth=1,
             label=f'Particle {particle_id}')
    
    # Add vertical line at beaching time if particle is beached
    beached_mask = particle_data['DistanceToShore'] < 0.001  # 1 meter threshold
    if beached_mask.any():
        beached_step = particle_data.loc[beached_mask.idxmax(), 'Step']
        ax2.axvline(x=beached_step, color=color, linestyle='--', alpha=0.3)

ax2.set_xlabel('Time Step')
ax2.set_ylabel('Distance to Shore (meters)')
ax2.set_title('Distance to Shore Over Time for All Particles\nVertical lines indicate beaching time')
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

# Calculate some statistics for each particle
stats = []
for particle_id in particle_ids:
    particle_data = df[df['ParticleID'] == particle_id]
    
    # Basic distance calculations
    total_distance = np.sum(np.sqrt(
        np.diff(particle_data['Longitude'])**2 + 
        np.diff(particle_data['Latitude'])**2
    ))
    
    # Velocity calculations
    velocities = np.sqrt(
        particle_data['VelocityU']**2 + 
        particle_data['VelocityV']**2
    )
    avg_velocity = np.mean(velocities)
    max_velocity = np.max(velocities)
    velocity_std = np.std(velocities)
    
    # Shore distance calculations
    min_shore_dist = particle_data['DistanceToShore'].min()
    avg_shore_dist = particle_data['DistanceToShore'].mean()
    shore_dist_std = particle_data['DistanceToShore'].std()
    
    # Time calculations
    total_time = particle_data['Step'].max() - particle_data['Step'].min()
    
    # Direction calculations
    dx = np.diff(particle_data['Longitude'])
    dy = np.diff(particle_data['Latitude'])
    angles = np.arctan2(dy, dx) * 180 / np.pi
    angle_std = np.std(angles)  # Standard deviation of direction changes
    
    # Check if particle is beached
    beached_mask = particle_data['DistanceToShore'] < 0.001  # 1 meter threshold
    beached_step = particle_data.loc[beached_mask.idxmax(), 'Step'] if beached_mask.any() else None
    
    # Calculate net displacement
    start_pos = np.array([particle_data.iloc[0]['Latitude'], particle_data.iloc[0]['Longitude']])
    end_pos = np.array([particle_data.iloc[-1]['Latitude'], particle_data.iloc[-1]['Longitude']])
    net_displacement = np.linalg.norm(end_pos - start_pos)
    
    # Calculate straightness index (net displacement / total distance)
    straightness = net_displacement / total_distance if total_distance > 0 else 0
    
    stats.append({
        'ParticleID': particle_id,
        'Total Distance': total_distance,
        'Net Displacement': net_displacement,
        'Straightness Index': straightness,  # 1 = straight line, 0 = very tortuous
        'Average Velocity': avg_velocity,
        'Max Velocity': max_velocity,
        'Velocity Std Dev': velocity_std,
        'Min Distance to Shore': min_shore_dist,
        'Avg Distance to Shore': avg_shore_dist,
        'Shore Distance Std Dev': shore_dist_std,
        'Total Time Steps': total_time,
        'Direction Change Std Dev': angle_std,  # Higher values indicate more erratic movement
        'Beached': beached_mask.any(),
        'Beaching Time Step': beached_step,
        'Start Latitude': particle_data.iloc[0]['Latitude'],
        'Start Longitude': particle_data.iloc[0]['Longitude'],
        'End Latitude': particle_data.iloc[-1]['Latitude'],
        'End Longitude': particle_data.iloc[-1]['Longitude']
    })

# Convert to DataFrame and save statistics
stats_df = pd.DataFrame(stats)

# Calculate some aggregate statistics
print("\nAggregate Statistics:")
print(f"Total number of particles: {len(particle_ids)}")
print(f"Total simulation steps: {df['Step'].max() + 1} steps")
print(f"Number of beached particles: {stats_df['Beached'].sum()}")
print(f"Percentage of beached particles: {(stats_df['Beached'].sum() / len(particle_ids) * 100):.2f}%")
print(f"Average straightness index: {stats_df['Straightness Index'].mean():.3f} (unitless)")
print(f"Average total distance: {stats_df['Total Distance'].mean():.3f} (degrees)")
print(f"Average net displacement: {stats_df['Net Displacement'].mean():.3f} (degrees)")
print(f"Average velocity: {stats_df['Average Velocity'].mean():.3f} m/s") # Assuming m/s based on context
print(f"Average time to beach (for beached particles): {stats_df[stats_df['Beached']]['Beaching Time Step'].mean():.1f} steps")

# Save detailed statistics
stats_df.to_csv('particle_statistics.csv', index=False)
print("\nDetailed statistics saved to 'particle_statistics.csv'") 