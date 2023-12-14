import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from collections import defaultdict

def initialize_boids(N):
    # Initialize Boids with random positions and zero velocities
    positions = np.random.uniform(-1, 1, size=(N, 3))
    velocities = np.zeros((N, 3))
    return positions, velocities

def spatial_hashing(positions, cell_size):
    # Perform spatial hashing to efficiently find neighbors
    hash_table = defaultdict(list)
    for i, pos in enumerate(positions):
        cell_key = tuple((pos / cell_size).astype(int))
        hash_table[cell_key].append(i)
    return hash_table

def get_neighbors_spatial_hashing(positions, i, hash_table, cell_size, do):
    # Get indices of Boids within the specified radius around the i-th Boid using spatial hashing
    cell_key = tuple((positions[i] / cell_size).astype(int))
    neighbors_indices = set()
    for x in range(cell_key[0] - 1, cell_key[0] + 2):
        for y in range(cell_key[1] - 1, cell_key[1] + 2):
            for z in range(cell_key[2] - 1, cell_key[2] + 2):
                cell = (x, y, z)
                neighbors_indices.update(hash_table[cell])
    neighbors_indices.discard(i)  # Remove the Boid itself
    return list(neighbors_indices)

def update_boids(positions, velocities, hash_table, cell_size, do, dc, l0, l1, l2, l3, l4, vmax):
    # Update velocities and positions of Boids based on the rules
    for i in range(len(positions)):
        w1 = rule1(positions, i, hash_table, cell_size, do)
        w2 = rule2(positions, velocities, i, hash_table, cell_size, do)
        w3 = rule3(positions, i, hash_table, cell_size, dc)
        w4 = rule4(positions[i], vmax)

        velocities[i] = l0 * velocities[i] + l1 * w1 + l2 * w2 + l3 * w3 + l4 * w4

        # Limit velocity to vmax
        norm_v = np.linalg.norm(velocities[i])
        if norm_v > vmax:
            velocities[i] = (vmax / norm_v) * velocities[i]

        # Update position
        positions[i] = positions[i] + velocities[i]

    return positions, velocities

def rule1(positions, i, hash_table, cell_size, do):
    # Rule 1: Move towards the center of mass of neighboring Boids using spatial hashing
    neighbors = get_neighbors_spatial_hashing(positions, i, hash_table, cell_size, do)
    if len(neighbors) == 0:
        return np.zeros(3)
    center_of_mass = np.mean(positions[neighbors], axis=0)
    return (center_of_mass - positions[i])

def rule2(positions, velocities, i, hash_table, cell_size, do):
    # Rule 2: Align velocity with the average velocity of neighboring Boids using spatial hashing
    neighbors = get_neighbors_spatial_hashing(positions, i, hash_table, cell_size, do)
    if len(neighbors) == 0:
        return np.zeros(3)
    avg_velocity = np.mean(velocities[neighbors], axis=0)
    return avg_velocity - velocities[i]

def rule3(positions, i, hash_table, cell_size, dc):
    # Rule 3: Avoid collisions with close neighbors using spatial hashing
    neighbors = get_neighbors_spatial_hashing(positions, i, hash_table, cell_size, dc)
    if len(neighbors) == 0:
        return np.zeros(3)
    avoidance_vector = np.mean(positions[neighbors] - positions[i], axis=0)
    return -avoidance_vector

def rule4(position, vmax):
    # Rule 4: Keep Boids within the cube
    return np.array([0 if abs(coord) <= 1 else np.sign(coord) * (abs(coord) - 1) for coord in position])

def visualize_boids(positions, ax):
    # Visualize Boids in 3D
    ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2], c='blue', marker='o')

# Parameters
N = 500
vmax = 0.03
do = 0.02
dc = 0.01
l0, l1, l2, l3, l4 = 0.31, 0.001, 1.2, 2, 0.01
cell_size = 0.2  # Adjust the cell size as needed

# Initialize Boids
positions, velocities = initialize_boids(N)

# Visualize initial state
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
visualize_boids(positions, ax)

# Simulation loop
hash_table = spatial_hashing(positions, cell_size)
for step in range(100):  # Adjust the number of steps as needed
    positions, velocities = update_boids(positions, velocities, hash_table, cell_size, do, dc, l0, l1, l2, l3, l4, vmax)
    ax.clear()
    visualize_boids(positions, ax)
    plt.pause(0.1)  # Adjust the pause time as needed

plt.show()
