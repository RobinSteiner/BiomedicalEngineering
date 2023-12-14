import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def initialize_boids(N):
    # Initialize Boids with random positions and zero velocities
    positions = np.random.uniform(-1, 1, size=(N, 3))
    velocities = np.zeros((N, 3))
    return positions, velocities

def update_boids(positions, velocities, do, dc, l0, l1, l2, l3, l4, vmax):
    # Update velocities and positions of Boids based on the rules
    for i in range(len(positions)):
        w1 = rule1(positions, i, do)
        w2 = rule2(positions, velocities, i, do)
        w3 = rule3(positions, i, dc)
        w4 = rule4(positions[i], vmax)

        velocities[i] = l0 * velocities[i] + l1 * w1 + l2 * w2 + l3 * w3 + l4 * w4

        # Limit velocity to vmax
        norm_v = np.linalg.norm(velocities[i])
        if norm_v > vmax:
            velocities[i] = (vmax / norm_v) * velocities[i]

        # Update position
        positions[i] = positions[i] + velocities[i]

    return positions, velocities

def rule1(positions, i, do):
    # Rule 1: Move towards the center of mass of neighboring Boids
    neighbors = get_neighbors(positions, i, do)
    if len(neighbors) == 0:
        return np.zeros(3)
    center_of_mass = np.mean(neighbors, axis=0)
    return (center_of_mass - positions[i])

def rule2(positions, velocities, i, do):
    # Rule 2: Align velocity with the average velocity of neighboring Boids
    neighbors = get_neighbors(positions, i, do)
    if len(neighbors) == 0:
        return np.zeros(3)
    avg_velocity = np.mean(velocities[neighbors], axis=0)
    return avg_velocity - velocities[i]

def rule3(positions, i, dc):
    # Rule 3: Avoid collisions with close neighbors
    neighbors = get_neighbors(positions, i, dc)
    if len(neighbors) == 0:
        return np.zeros(3)
    avoidance_vector = np.mean(positions[neighbors] - positions[i], axis=0)
    return -avoidance_vector

def rule4(position, vmax):
    # Rule 4: Keep Boids within the cube
    return np.array([0 if abs(coord) <= 1 else np.sign(coord) * (abs(coord) - 1) for coord in position])

def get_neighbors(positions, i, radius):
    # Get indices of Boids within the specified radius around the i-th Boid
    distances = np.linalg.norm(positions - positions[i], axis=1)
    neighbors_indices = np.where((distances > 0) & (distances <= radius))[0]
    return neighbors_indices

def visualize_boids(positions, ax):
    # Visualize Boids in 3D
    ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2], c='blue', marker='o')

# Parameters
N = 500
vmax = 0.03
do = 0.02
dc = 0.01
l0, l1, l2, l3, l4 = 0.31, 0.001, 1.2, 2, 0.01

# Initialize Boids
positions, velocities = initialize_boids(N)

# Visualize initial state
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
visualize_boids(positions, ax)

# Simulation loop
for step in range(1000):  # Adjust the number of steps as needed
    positions, velocities = update_boids(positions, velocities, do, dc, l0, l1, l2, l3, l4, vmax)
    ax.clear()
    visualize_boids(positions, ax)
    plt.pause(0.1)  # Adjust the pause time as needed

plt.show()
