import h5py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# === Load shoreline mask ===
with h5py.File("../../../data/shore_distance.h5", "r") as f:
    shore_lat = f["lat"][:]
    shore_lon = f["lon"][:]
    shore_mask = f["mask"][:].reshape(len(shore_lat), len(shore_lon))
    shore_distance = f["distance"][:].reshape(len(shore_lat), len(shore_lon))

shore_lon_grid, shore_lat_grid = np.meshgrid(shore_lon, shore_lat)

# === Load velocity grid ===
with h5py.File("../../../data/grid.h5", "r") as f:
    grid_lat = f["latitude"][:]
    grid_lon = f["longitude"][:]

grid_lon_grid, grid_lat_grid = np.meshgrid(grid_lon, grid_lat)

# === Define input files and styles ===
files = [
    ("snap_trajectory_dualcondition.csv", "Snap + DualCondition", "red"),
    ("freeslip_trajectory_distancebased.csv", "FreeSlip + DistanceBased", "blue"),
    ("snap_trajectory_velocitybased.csv", "Snap + VelocityBased", "green"),
    ("partial_trajectory_dualcondition.csv", "PartialSlip + DualCondition", "orange")
]

# === Create figure with 2 subplots ===
fig = plt.figure(figsize=(15, 10))
ax1 = plt.subplot(121, projection=ccrs.PlateCarree())
ax2 = plt.subplot(122)

# Determine map zoom window
all_lons, all_lats, all_beach_steps = [], [], []

for filename, _, _ in files:
    df = pd.read_csv(filename)
    all_lons.extend(df['Longitude'])
    all_lats.extend(df['Latitude'])

    if df['Beached'].eq('Yes').any():
        step = df[df['Beached'] == 'Yes']['Step'].iloc[0]
        all_beach_steps.append(step)

# Map extent
lon_min, lon_max = min(all_lons), max(all_lons)
lat_min, lat_max = min(all_lats), max(all_lats)
lon_pad = (lon_max - lon_min) * 0.1
lat_pad = (lat_max - lat_min) * 0.1
ax1.set_extent([lon_min - lon_pad, lon_max + lon_pad, lat_min - lat_pad, lat_max + lat_pad])

# Crop shore_distance grid
bbox = [lon_min - lon_pad, lon_max + lon_pad, lat_min - lat_pad, lat_max + lat_pad]
lon_mask = (shore_lon >= bbox[0]) & (shore_lon <= bbox[1])
lat_mask = (shore_lat >= bbox[2]) & (shore_lat <= bbox[3])
shore_lon_crop = shore_lon[lon_mask]
shore_lat_crop = shore_lat[lat_mask]
shore_distance_crop = shore_distance[np.ix_(lat_mask, lon_mask)]
shore_lon_grid_crop, shore_lat_grid_crop = np.meshgrid(shore_lon_crop, shore_lat_crop)

# Plot shore grid
c = ax1.pcolormesh(shore_lon_grid_crop, shore_lat_grid_crop, shore_distance_crop,
                   cmap='Blues_r', shading='auto', transform=ccrs.PlateCarree())
plt.colorbar(c, ax=ax1, label='Distance to Shore (meters)')

# Plot shoreline mask
ax1.contourf(shore_lon_grid, shore_lat_grid, shore_mask, levels=[0.5, 1], colors='lightgray', alpha=0.8, transform=ccrs.PlateCarree())

# Plot velocity grid
ax1.plot(grid_lon_grid, grid_lat_grid, color='lightblue', linewidth=0.8, transform=ccrs.PlateCarree())
ax1.plot(grid_lon_grid.T, grid_lat_grid.T, color='lightblue', linewidth=0.8, transform=ccrs.PlateCarree())

# === Plot each trajectory on map and time-series ===
for filename, label, color in files:
    df = pd.read_csv(filename)

    # Map trajectory
    ax1.plot(df['Longitude'], df['Latitude'], '-', label=label, color=color, transform=ccrs.PlateCarree())
    ax1.scatter(df['Longitude'].iloc[0], df['Latitude'].iloc[0], color=color, s=50, marker='o', transform=ccrs.PlateCarree())

    # Distance-to-shore time series
    ax2.plot(df['Step'], df['DistanceToShore'], label=label, color=color)

    if df['Beached'].eq('Yes').any():
        beach_idx = df[df['Beached'] == 'Yes'].index[0]
        beach_step = df.loc[beach_idx, 'Step']
        beach_dist = df.loc[beach_idx, 'DistanceToShore']

        ax2.axvline(x=beach_step, linestyle='--', color=color, alpha=0.5)
        ax2.text(beach_step, beach_dist + 2,
                 f"{label}\nBeached at {beach_dist:.1f} m",
                 rotation=90, color=color,
                 verticalalignment='bottom', horizontalalignment='center')

# Zoom in around beaching steps
if all_beach_steps:
    min_step = max(min(all_beach_steps) - 10, 0)
    max_step = max(all_beach_steps) + 20
    ax2.set_xlim(min_step, max_step)

# Final styling
ax1.set_title('Particle Trajectories Near Coastline')
ax1.legend()

ax2.set_title('Distance to Shore Over Time (Zoomed)')
ax2.set_xlabel('Time Step')
ax2.set_ylabel('Distance to Shore (meters)')
ax2.grid(True)
# === Build dynamic legend with beaching info ===
legend_labels = []

for filename, label, color in files:
    df = pd.read_csv(filename)

    if df['Beached'].eq('Yes').any():
        beach_idx = df[df['Beached'] == 'Yes'].index[0]
        beach_step = df.loc[beach_idx, 'Step']
        beach_dist = df.loc[beach_idx, 'DistanceToShore']
        legend_labels.append((f"{label} (Step {beach_step}, {beach_dist:.1f} m)", color))
    else:
        legend_labels.append((f"{label} (No beaching)", color))

# Add updated legend with full info
handles = [plt.Line2D([0], [0], color=color, lw=2) for _, color in legend_labels]
labels = [text for text, _ in legend_labels]
ax2.legend(handles, labels, loc='best')

plt.tight_layout()
plt.savefig("multi_trajectory_zoom_annotated.png", dpi=300)
plt.show()
