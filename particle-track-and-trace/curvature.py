import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load your CSV
df = pd.read_csv('single_particle_trajectory.csv')

# Extract x (lon), y (lat)
x = df['Longitude'].values
y = df['Latitude'].values

# Time step (you know dt = 3600 sec in your sim)
dt = 3600.0  # seconds per step, adjust if needed

# First derivative: velocity (numerical diff)
dx = np.gradient(x, dt)
dy = np.gradient(y, dt)

# Second derivative: acceleration
ddx = np.gradient(dx, dt)
ddy = np.gradient(dy, dt)

# Curvature
curvature = (dx * ddy - dy * ddx) / ( (dx**2 + dy**2)**1.5 )

# Optional: clean up NaNs / infs at start and end
curvature = np.nan_to_num(curvature)

# Plot curvature over time step
plt.figure(figsize=(12, 6))
plt.plot(df['Step'], curvature, label='Curvature κ(t)')
plt.xlabel('Time Step')
plt.ylabel('Curvature κ(t)')
plt.title('Path Curvature Over Time')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('curvature_over_time.png', dpi=300)
plt.show()

# Compute summary statistics
mean_curvature = np.mean(np.abs(curvature))
max_curvature = np.max(np.abs(curvature))
curvature_variance = np.var(curvature)

print(f"Mean curvature: {mean_curvature:.6f}")
print(f"Max curvature: {max_curvature:.6f}")
print(f"Curvature variance: {curvature_variance:.6f}")