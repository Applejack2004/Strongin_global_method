import matplotlib.pyplot as plt
import numpy as np
import os
import sys

silent_mode = False
task_index = ""
if len(sys.argv) >= 3 and sys.argv[1] == "silent":
    silent_mode = True
    task_index = sys.argv[2]

data_file_path = "function_plot_data_2d.txt"
points_file_path = "trial_points.txt"

def main():
    if not os.path.exists(data_file_path) or not os.path.exists(points_file_path):
        return

    with open(data_file_path, 'r') as f:
        header = f.readline().split()
        num_constrs = int(header[0])
        opt_y1, opt_y2 = float(header[1]), float(header[2])

    grid_data = np.loadtxt(data_file_path, skiprows=1)
    y1_vals = np.unique(grid_data[:, 0])
    y2_vals = np.unique(grid_data[:, 1])
    Y1, Y2 = np.meshgrid(y1_vals, y2_vals)
    
    F = grid_data[:, 2].reshape(len(y1_vals), len(y2_vals)).T
    G =[]
    for i in range(num_constrs):
        G.append(grid_data[:, 3+i].reshape(len(y1_vals), len(y2_vals)).T)

    trials = np.loadtxt(points_file_path)
    if trials.ndim == 1: trials = trials.reshape(1, -1)
    
    ty1 = trials[:, 0] * 4.0
    ty2 = -1.0 + trials[:, 1] * 4.0
    tv = trials[:, 3].astype(int)

    plt.figure(figsize=(10, 9))
    levels = np.linspace(np.min(F), np.max(F), 40)
    plt.contour(Y1, Y2, F, levels=levels, cmap='viridis', alpha=0.3, linewidths=0.5)

    colors =['red', 'green', 'magenta']
    for i in range(num_constrs):
        plt.contour(Y1, Y2, G[i], levels=[0], colors=colors[i], linewidths=2.5)

    infeasible_mask = np.zeros_like(F, dtype=bool)
    for i in range(num_constrs):
        infeasible_mask |= (G[i] > 0)
    
    plt.imshow(infeasible_mask, extent=(y1_vals.min(), y1_vals.max(), y2_vals.min(), y2_vals.max()), origin='lower', cmap='Reds', alpha=0.1, aspect='auto')

    mask_feasible = (tv == num_constrs + 1)
    plt.scatter(ty1[~mask_feasible], ty2[~mask_feasible], color='black', s=2, alpha=0.4, label='Infeasible Trials')
    plt.scatter(ty1[mask_feasible], ty2[mask_feasible], color='blue', s=12, alpha=0.8, edgecolors='white', linewidth=0.3, label='Feasible Trials')
    plt.scatter(opt_y1, opt_y2, color='gold', marker='*', s=300, edgecolors='black', zorder=10, label='Found Minimum')

    plt.title("Multidimensional Constrained Search Map\nContour Lines and Constraint Boundaries")
    plt.xlabel("y1")
    plt.ylabel("y2")
    plt.legend(loc='upper right')
    plt.grid(True, linestyle=':', alpha=0.5)
    plt.xlim(0, 4)
    plt.ylim(-1, 3)
    plt.tight_layout()

    if silent_mode:
        os.makedirs("plots_output", exist_ok=True)
        plt.savefig(f"plots_output/2d_constr_task_{task_index}.png", dpi=150)
    else:
        plt.show()

if __name__ == "__main__":
    main()