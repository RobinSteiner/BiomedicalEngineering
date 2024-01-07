import numpy as np
from collections import defaultdict
from vispy import app
from vispy import scene
from vispy.scene import visuals


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


def get_neighbors(positions, i, hash_table, cell_size, d):
    # Get indices of Boids within the specified radius around the i-th Boid using spatial hashing
    cell_key = tuple((positions[i] / cell_size).astype(int))
    neighbors_indices = set()
    for x in range(cell_key[0] - 1, cell_key[0] + 1):
        for y in range(cell_key[1] - 1, cell_key[1] + 1):
            for z in range(cell_key[2] - 1, cell_key[2] + 1):
                cell = (x, y, z)
                if hash_table[cell] != i and np.linalg.norm(positions[i] - positions[hash_table[cell]]) <= d:
                    neighbors_indices.update(hash_table[cell])
    return list(neighbors_indices)


def update_boids(positions, velocities, hash_table, cell_size, do, dc, l0, l1, l2, l3, l4, vmax):
    # Update velocities and positions of Boids based on the rules
    for i in range(len(positions)):
        neighbors = get_neighbors(positions, i, hash_table, cell_size, do)
        neighbors1 = get_neighbors(positions, i, hash_table, cell_size, dc)

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
    # Rule 1: Move towards the center of mass of neighboring Boids
    if len(neighbors) == 0:
        return np.zeros(3)
    center_of_mass = np.mean(positions[neighbors], axis=0)
    return center_of_mass - positions[i]


def rule2(velocities, neighbors, i):
    # Rule 2: Align velocity with the average velocity of neighboring Boids
    if len(neighbors) == 0:
        return np.zeros(3)
    avg_velocity = np.mean(velocities[neighbors], axis=0)
    return avg_velocity


def rule3(positions, neighbors, i):
    # Rule 3: Avoid collisions with close neighbors
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


class BoidsVisual(scene.SceneCanvas):
    def __init__(self, positions, velocities):
        scene.SceneCanvas.__init__(self, keys='interactive', size=(800, 600), bgcolor='white')
        self.unfreeze()

        self.positions = positions
        self.velocities = velocities

        self.lines = visuals.Line(width=2)
        self.heads = visuals.Markers()
        self.update_boids()

        self.view = self.central_widget.add_view()
        self.view.add(self.lines)
        self.view.add(self.heads)
        self.view.camera = scene.cameras.TurntableCamera(up='z')

        self.timer = app.Timer('auto', self.on_timer)
        self.timer.start()

        self.show()

    def update_boids(self):
        # Calculate color of Boid
        min_pos = self.positions.min(axis=0)
        max_pos = self.positions.max(axis=0)
        norm_positions = (self.positions - min_pos) / (max_pos - min_pos)
        colors = np.hstack((norm_positions, np.ones((len(norm_positions), 1))))
        segment_colors = np.empty((len(self.positions) * 2, 4), dtype=np.float32)
        segment_colors[0::2] = colors
        segment_colors[1::2] = colors

        # Separate positions and endpoints into two (N, 3) arrays
        pos_data = np.zeros((len(self.positions) * 2, 3), dtype=np.float32)
        pos_data[0::2] = self.positions
        pos_data[1::2] = self.positions - self.velocities

        # Update the visual
        self.heads.set_data(self.positions, edge_color=None, face_color=colors, size=4)
        self.lines.set_data(pos=pos_data, connect='segments', color=segment_colors)

    def on_timer(self, event):
        global positions, velocities, hash_table, cell_size, do, dc, l0, l1, l2, l3, l4, v_max
        hash_table = spatial_hashing(positions, cell_size)
        positions, velocities = update_boids(positions, velocities, hash_table, cell_size, do, dc, l0, l1, l2, l3, l4,
                                             v_max)

        self.update_boids()
        self.update()


# Parameters
N = 500
v_max = 0.03
do = 0.2
dc = 0.1
l0, l1, l2, l3, l4 = 0.31, 0.1, 10, 2, 2
cell_size = 0.1

# Initialize Boids
positions, velocities = initialize_boids(N)
hash_table = spatial_hashing(positions, cell_size)

c = BoidsVisual(positions, velocities)
app.run()
