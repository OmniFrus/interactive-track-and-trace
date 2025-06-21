import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load your CSV
df = pd.read_csv("all_particles_trajectory.csv")

# Select one particle for now (ParticleID = 0)
particle_id = 100
df_particle = df[df["ParticleID"] == particle_id].copy()

# Time step in seconds (your simulation dt = 3600 seconds)
dt = 3600

# Extract coordinates
x = df_particle["Longitude"].values
y = df_particle["Latitude"].values

# Compute velocity (Jacobian)
vx = np.gradient(x, dt)
vy = np.gradient(y, dt)

# Compute acceleration (Hessian)
ax = np.gradient(vx, dt)
ay = np.gradient(vy, dt)

# Compute norms
v_norm = np.sqrt(vx**2 + vy**2)
a_norm = np.sqrt(ax**2 + ay**2)

# Global metrics
mean_velocity = np.mean(v_norm)
mean_acceleration = np.mean(a_norm)

print(f"Particle {particle_id} — Mean velocity (deg/s): {mean_velocity:.6f}")
print(f"Particle {particle_id} — Mean acceleration (deg/s²): {mean_acceleration:.6f}")

# Plot trajectory with velocity and acceleration
plt.figure(figsize=(12, 6))

plt.subplot(1, 2, 1)
plt.plot(x, y, 'b-', label='Trajectory')
plt.quiver(x, y, vx, vy, color='r', scale=0.05, label='Velocity')
plt.title("Trajectory and Jacobian (velocity)")
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(df_particle["Step"], v_norm, label='Velocity norm')
plt.plot(df_particle["Step"], a_norm, label='Acceleration norm')
plt.xlabel("Step")
plt.ylabel("Norm")
plt.title("Jacobian and Hessian norms over time")
plt.legend()

plt.tight_layout()
plt.savefig("jacobian_hessian_particle.png", dpi=300)
plt.show()