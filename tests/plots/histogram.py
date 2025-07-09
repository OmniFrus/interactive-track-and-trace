import pandas as pd
import matplotlib.pyplot as plt

# === File paths ===
velocity_file = 'snap_velocity_trajectories.csv'
dual_file = 'snap_dualcondition_trajectories.csv'

# === Load and filter each ===
def load_beached_distances(filepath):
    df = pd.read_csv(filepath)
    if 'Beached' in df.columns:
        if df['Beached'].dtype == object:
            df = df[df['Beached'].str.lower() == 'yes']
        else:
            df = df[df['Beached'] == True]
    return df['DistanceToShore'].dropna()

# Get data
velocity_distances = load_beached_distances(velocity_file)
dual_distances = load_beached_distances(dual_file)

# === Plot combined histogram ===
plt.figure(figsize=(10, 6))
bins = range(0, 22000, 500)

plt.hist(velocity_distances, bins=bins, alpha=0.6, color='skyblue', label='Snap + Velocity-Based')
plt.hist(dual_distances, bins=bins, alpha=0.6, color='orange', label='Snap + Dual-Condition')

plt.xlabel('Distance to Shore at Beaching (meters)')
plt.ylabel('Number of Particles')
plt.title('Beaching Distance Distribution: Snap Configurations')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('snap_config_beaching_histogram.png', dpi=300)
plt.show()