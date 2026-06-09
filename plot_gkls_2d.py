import matplotlib.pyplot as plt
import numpy as np
import os
import sys

silent_mode = False
task_index = ""
if len(sys.argv) >= 3 and sys.argv[1] == "silent":
    silent_mode = True
    task_index = sys.argv[2]

def main():
    data_path = "function_plot_data_2d.txt"
    points_path = "trial_points.txt"
    if not os.path.exists(data_path): return

    with open(data_path, 'r') as f:
        header = f.readline().split()
        num_constrs = int(header[0])
        opt_x, opt_y = float(header[1]), float(header[2])

    grid_data = np.loadtxt(data_path, skiprows=1)
    y1_vals = np.unique(grid_data[:, 0])
    y2_vals = np.unique(grid_data[:, 1])
    Y1, Y2 = np.meshgrid(y1_vals, y2_vals)
    F = grid_data[:, 2].reshape(len(y1_vals), len(y2_vals)).T

    plt.figure(figsize=(10, 8))
    contour = plt.contourf(Y1, Y2, F, levels=30, cmap='viridis', alpha=0.7)
    plt.colorbar(contour, label='f(x)')

    for i in range(num_constrs):
        G = grid_data[:, 3+i].reshape(len(y1_vals), len(y2_vals)).T
        plt.contour(Y1, Y2, G, levels=[0], colors='red', linewidths=2)
        plt.contourf(Y1, Y2, G, levels=[0, 1e9], colors='red', alpha=0.1)

    if os.path.exists(points_path):
        pts = np.loadtxt(points_path)
        if pts.ndim == 1: pts = pts.reshape(1, -1)
        tx = -1.0 + 2.0 * pts[:, 0]
        ty = -1.0 + 2.0 * pts[:, 1]
        tv = pts[:, 3].astype(int)
        mask = (tv == num_constrs + 1)
        plt.scatter(tx[~mask], ty[~mask], c='black', s=5, alpha=0.3, label='Infeasible')
        plt.scatter(tx[mask], ty[mask], c='blue', s=10, alpha=0.6, label='Feasible')

    plt.scatter(opt_x, opt_y, color='gold', marker='*', s=250, edgecolors='black', label='Global Opt')
    plt.title(f"GKLS Constrained Landscape (Task {task_index})" if silent_mode else "GKLS Constrained Landscape")
    plt.legend()
    plt.xlim(-1, 1); plt.ylim(-1, 1)
    plt.tight_layout()

    if silent_mode:
        os.makedirs("plots_output", exist_ok=True)
        plt.savefig(f"plots_output/gkls_2d_task_{task_index}.png", dpi=150)
    else:
        plt.show()

if __name__ == "__main__":
    main()