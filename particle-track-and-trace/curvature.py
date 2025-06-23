import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# === LOAD CSV ===
df = pd.read_csv("single_particle_trajectory.csv")

# === PREPARE DATA ===
# 1 step = 1 hour
df['time_hours'] = df['Step'] * 1.0

# x = Longitude, y = Latitude
x = df['Longitude'].values
y = df['Latitude'].values
t = df['time_hours'].values

# === COMPUTE DERIVATIVES ===
dx_dt = np.gradient(x, t)
dy_dt = np.gradient(y, t)
d2x_dt2 = np.gradient(dx_dt, t)
d2y_dt2 = np.gradient(dy_dt, t)

# === COMPUTE CURVATURE ===
epsilon = 1e-10
curvature = np.abs(dx_dt * d2y_dt2 - dy_dt * d2x_dt2) / ((dx_dt**2 + dy_dt**2)**1.5 + epsilon)

# === FFT ===
fft_result = np.fft.fft(curvature)
fft_freqs = np.fft.fftfreq(len(curvature), d=(t[1]-t[0]))  # d = 1 hour
fft_magnitude = np.abs(fft_result)

# === PLOT CURVATURE vs TIME ===
plt.figure(figsize=(12,5))
plt.plot(t, curvature)
plt.title("Curvature vs Time")
plt.xlabel("Time (hours)")
plt.ylabel("Curvature κ(t)")
plt.grid()
plt.show()

# === PLOT FFT ===
positive_freqs = fft_freqs[fft_freqs >= 0]
positive_magnitude = fft_magnitude[fft_freqs >= 0]

# Only positive frequencies > 0
nonzero_freqs = positive_freqs[positive_freqs > 0]
nonzero_magnitude = positive_magnitude[positive_freqs > 0]

# Only keep frequencies > threshold
freq_threshold = 0.01  # (1/hour), for example: skip ultra-low frequencies

valid_mask = (nonzero_freqs > freq_threshold)
valid_freqs = nonzero_freqs[valid_mask]
valid_magnitude = nonzero_magnitude[valid_mask]

dominant_index = np.argmax(valid_magnitude)
dominant_freq = valid_freqs[dominant_index]
dominant_value = valid_freqs[dominant_index]

# Plot
plt.figure(figsize=(12,6))
plt.plot(valid_freqs, valid_magnitude, marker='o')
plt.title("Fourier Spectrum of Path Curvature (no DC)")
plt.xlabel("Frequency (1/hour)")
plt.ylabel("Magnitude")
plt.grid()
plt.annotate(f'Dominant freq = {dominant_freq:.4f} (1/hour)', xy=(dominant_freq, dominant_value),
             xytext=(dominant_freq, dominant_value*1.1),
             arrowprops=dict(facecolor='red', shrink=0.05))
plt.show()

# === PRINT DOMINANT FREQUENCY ===
print(f"Dominant frequency: {dominant_freq:.4f} (1/hour)")