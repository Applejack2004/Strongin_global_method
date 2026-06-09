import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from scipy.interpolate import griddata

POINTS_FILE = 'grishagin_points.txt'
LEVELS = 30  
GRID_RESOLUTION = 300 
OUTPUT_IMAGE_FILE = 'grishagin_visualization.png' 

silent_mode = False
task_index = ""
if len(sys.argv) >= 5 and sys.argv[3] == "silent":
    silent_mode = True
    task_index = sys.argv[4]

def load_points(filename):
    points_data =[]
    with open(filename, 'r') as f:
        for line in f:
            try:
                parts = line.strip().replace(',', ' ').split()
                if len(parts) < 3:
                    continue
                points_data.append(list(map(float, parts[:3])))
            except (ValueError, IndexError):
                continue

    if not points_data:
        return np.array([]), np.array([]), np.array([])

    points_data = np.array(points_data)
    return points_data[:, 0], points_data[:, 1], points_data[:, 2]

if __name__ == '__main__':
    if not os.path.exists(POINTS_FILE):
        exit(1)

    actual_min_coords = None
    if len(sys.argv) >= 3:
        try:
            actual_x = float(sys.argv[1])
            actual_y = float(sys.argv[2])
            actual_min_coords = (actual_x, actual_y)
        except ValueError:
            pass

    xs, ys, zs = load_points(POINTS_FILE)

    def filter_points(xs, ys, zs):
        mask = (xs >= 0) & (xs <= 1) & (ys >= 0) & (ys <= 1)
        return xs[mask], ys[mask], zs[mask]

    xs, ys, zs = filter_points(xs, ys, zs)

    if len(xs) < 4:
        exit(0)

    grid_x, grid_y = np.meshgrid(
        np.linspace(0, 1, GRID_RESOLUTION),
        np.linspace(0, 1, GRID_RESOLUTION)
    )

    grid_z = griddata((xs, ys), zs, (grid_x, grid_y), method='cubic')

    plt.figure(figsize=(10, 10))  
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')

    contour = plt.contourf(grid_x, grid_y, grid_z, levels=LEVELS, cmap='viridis_r', alpha=0.8)
    plt.colorbar(contour, label='f(y1, y2) - interpolated value')
    plt.contour(grid_x, grid_y, grid_z, levels=contour.levels, colors='black', linewidths=0.5, alpha=0.5)

    plt.scatter(xs, ys, s=40, c='red', edgecolors='black', linewidth=0.5, label=f'Trial points ({len(xs)})', zorder=3)

    min_idx = np.argmin(zs)
    plt.scatter(xs[min_idx], ys[min_idx], s=200, c='cyan', marker='*', edgecolors='black', label=f'Found minimum: z={zs[min_idx]:.4f}', zorder=4)

    if actual_min_coords:
        plt.scatter(actual_min_coords[0], actual_min_coords[1], s=250, c='lime', marker='X', edgecolors='black', label='Actual minimum', zorder=5) 

    plt.xlabel('y1')
    plt.ylabel('y2')
    plt.title(f'Contour lines and search trajectory (Grishagin task {task_index})' if silent_mode else 'Contour lines and search trajectory (Grishagin)')
    plt.legend(loc='upper left')
    plt.axis('scaled')
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()

    if silent_mode:
        os.makedirs("plots_output", exist_ok=True)
        plt.savefig(f"plots_output/grishagin_task_{task_index}.png", dpi=150)
    else:
        plt.savefig(OUTPUT_IMAGE_FILE, dpi=150)
        plt.show()