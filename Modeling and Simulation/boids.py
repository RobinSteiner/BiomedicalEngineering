import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

def initialize_boids(N):
    # Initialize Boids with random positions and zero velocities
    positions = np.random.uniform(-1, 1, size=(N, 3))
    velocities = np.zeros((N, 3))
    return positions, velocities


def get_neighbors(positions, i, do):
    # Get indices of Boids within the specified radius around the i-th Boid using spatial hashing
    neighbors_indices = list()
    for x in range(len(positions)):
        if x != i and np.linalg.norm(positions[i] - positions[x]) <= do:
            neighbors_indices.append(x)
    return neighbors_indices


def update_boids(positions, velocities, do, dc, l0, l1, l2, l3, l4, vmax):
    # Update velocities and positions of Boids based on the rules
    for i in range(len(positions)):
        neighbors = get_neighbors(positions, i, do)
        neighbors1 = get_neighbors(positions, i, dc)

        w1 = rule1(positions, neighbors, i)
        w2 = rule2(velocities, neighbors, i)
        w3 = rule3(positions, neighbors1, i)
        w4 = rule4(positions, vmax, i)

        velocities[i] = l0 * velocities[i] + l1 * w1 + l2 * w2 + l3 * w3 + l4 * w4

        # Limit velocity to vmax
        norm_v = np.linalg.norm(velocities[i])
        if norm_v > vmax:
            velocities[i] = (vmax / norm_v) * velocities[i]

        # Update position
        positions[i] = positions[i] + velocities[i]

    return positions, velocities


def rule1(positions, neighbors, i):
    # Rule 1: Move towards the center of mass of neighboring Boids using spatial hashing
    if len(neighbors) == 0:
        return np.zeros(3)
    center_of_mass = np.mean(positions[neighbors], axis=0)
    return center_of_mass - positions[i]


def rule2(velocities, neighbors, i):
    # Rule 2: Align velocity with the average velocity of neighboring Boids using spatial hashing
    if len(neighbors) == 0:
        return np.zeros(3)
    avg_velocity = np.mean(velocities[neighbors], axis=0)
    return avg_velocity


def rule3(positions, neighbors, i):
    # Rule 3: Avoid collisions with close neighbors using spatial hashing
    if len(neighbors) == 0:
        return np.zeros(3)
    avoidance_vector = np.mean(positions[neighbors] - positions[i], axis=0)
    return -avoidance_vector


def rule4(positions, vmax, i):
    # Rule 4: Keep Boids within the cube
    return np.where(abs(positions[i]) > 1, -positions[i] / np.linalg.norm(positions[i]) * vmax, np.zeros(3))

def visualize_boids(positions, ax):
    # Visualize Boids in 3D
    ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2], c='blue', marker='o')


# Parameters
N = 300
v_max = 0.03
do = 0.2
dc = 0.1
l0, l1, l2, l3, l4 = 0.31, 0.001, 1.2, 2, 1

# Initialize Boids
positions, velocities = initialize_boids(N)

# Visualize initial state
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim([-1, 1])
ax.set_ylim([-1, 1])
ax.set_zlim([-1, 1])

for step in range(1000):  # Adjust the number of steps as needed
    positions, velocities = update_boids(positions, velocities, do, dc, l0, l1, l2, l3, l4, v_max)

    ax.clear()
    visualize_boids(positions, ax)

    # Keep the axes limits constant
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])

    plt.pause(0.0001)

plt.show()
