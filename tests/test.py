import numpy as np
import matplotlib.pyplot as plt

def get_shore_distance(lat, lon):
    # Flat coast at x=0 (longitude)
    return abs(lon)

def grad_shore_distance(lat, lon, delta=1e-4):
    # Numerical gradient as in your C++ code
    grad_lat = (get_shore_distance(lat + delta, lon) - get_shore_distance(lat - delta, lon)) / (2 * delta)
    grad_lon = (get_shore_distance(lat, lon + delta) - get_shore_distance(lat, lon - delta)) / (2 * delta)
    return grad_lat, grad_lon

def metre_to_degrees(m):
    # Dummy conversion for test (use your real function if needed)
    return m * 1e-5

def advect_freeslip(time, lat, lon, dt, u, v):
    shore_dist = get_shore_distance(lat, lon)
    shoreThreshold = 0.02  # meters
    epsilon = 1e-5
    # Use the same logic as your C++ code
    if shore_dist < shoreThreshold:
        grad_lat, grad_lon = grad_shore_distance(lat, lon)
        norm = np.sqrt(grad_lat**2 + grad_lon**2)
        if norm >= 1e-8:
            normal_x = grad_lon / norm
            normal_y = grad_lat / norm
            tangent_x = -normal_y
            tangent_y = normal_x
            v_t = u * tangent_x + v * tangent_y
            slip_u = v_t * tangent_x
            slip_v = v_t * tangent_y
            new_lon = lon + metre_to_degrees(slip_u * dt)
            new_lat = lat + metre_to_degrees(slip_v * dt)
            return new_lat, new_lon
        else:
            return lat, lon
    # If not near shore, use normal advection
    new_lon = lon + metre_to_degrees(u * dt)
    new_lat = lat + metre_to_degrees(v * dt)
    return new_lat, new_lon

def advect_partialslip(time, lat, lon, dt, u, v):
    shore_dist = get_shore_distance(lat, lon)
    shoreThreshold = 0.02  # meters
    epsilon = 1e-5
    slipRatio = 0.5
    if shore_dist < shoreThreshold:
        grad_lat, grad_lon = grad_shore_distance(lat, lon)
        norm = np.sqrt(grad_lat**2 + grad_lon**2)
        if norm >= 1e-8:
            normal_x = grad_lon / norm
            normal_y = grad_lat / norm
            tangent_x = -normal_y
            tangent_y = normal_x
            v_n = u * normal_x + v * normal_y
            v_t = u * tangent_x + v * tangent_y
            slip_u = (1.0 - slipRatio) * v_n * normal_x + v_t * tangent_x
            slip_v = (1.0 - slipRatio) * v_n * normal_y + v_t * tangent_y
            new_lon = lon + metre_to_degrees(slip_u * dt)
            new_lat = lat + metre_to_degrees(slip_v * dt)
            return new_lat, new_lon
        else:
            return lat, lon
    new_lon = lon + metre_to_degrees(u * dt)
    new_lat = lat + metre_to_degrees(v * dt)
    return new_lat, new_lon

# Test parameters
dt = 1
n_steps = 1000
lat, lon = 0.0, -0.1
u, v = 100, 100

positions_freeslip = []
positions_partialslip = []

# FreeSlip
lat_fs, lon_fs = lat, lon
for _ in range(n_steps):
    lat_fs, lon_fs = advect_freeslip(0, lat_fs, lon_fs, dt, u, v)
    positions_freeslip.append((lon_fs, lat_fs))

# PartialSlip
lat_ps, lon_ps = lat, lon
for _ in range(n_steps):
    lat_ps, lon_ps = advect_partialslip(0, lat_ps, lon_ps, dt, u, v)
    positions_partialslip.append((lon_ps, lat_ps))

positions_freeslip = np.array(positions_freeslip)
positions_partialslip = np.array(positions_partialslip)
plt.plot(positions_freeslip[:,0], positions_freeslip[:,1], label='FreeSlip')
plt.plot(positions_partialslip[:,0], positions_partialslip[:,1], label='PartialSlip')
plt.axvline(0, color='k', linestyle='--', label='Coastline')
plt.legend()
plt.xlabel('Longitude (x)')
plt.ylabel('Latitude (y)')
plt.title('Test: Particle Approaching a Flat Coast (C++-like logic)')
plt.show()