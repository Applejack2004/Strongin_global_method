import matplotlib.pyplot as plt
import numpy as np
import os
import sys

silent_mode = False
task_index = ""
if len(sys.argv) >= 3 and sys.argv[1] == "silent":
    silent_mode = True
    task_index = sys.argv[2]

data_file_path = "function_plot_data.txt"
points_file_path = "trial_points.txt"

def main():
    if not os.path.exists(data_file_path): return

    try:
        with open(data_file_path, 'r') as f:
            header = f.readline().split()
            if len(header) < 3: return
            num_constraints = int(header[0])
            found_x = float(header[1])
            found_y = float(header[2])
            data = np.loadtxt(f)
    except Exception:
        return

    if data.size == 0: return

    x = data[:, 0]
    y_obj = data[:, 1]
    constraints =[]
    for i in range(num_constraints):
        constraints.append(data[:, 2 + i])

    trial_x, trial_z, trial_v = [], [],[]
    has_trials = False
    if os.path.exists(points_file_path):
        try:
            p_data = np.loadtxt(points_file_path)
            if p_data.size > 0:
                if p_data.ndim == 1: p_data = p_data.reshape(1, -1)
                trial_x = p_data[:, 0]
                trial_z = p_data[:, 1]
                trial_v = p_data[:, 2].astype(int)
                has_trials = True
        except Exception:
            pass

    infeasible_mask = np.zeros_like(x, dtype=bool)
    if num_constraints > 0:
        for g_vals in constraints:
            infeasible_mask = infeasible_mask | (g_vals > 0)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)

    ax1.set_title(f"Full Problem View: Objective, Constraints ({num_constraints}) and All Trials")
    ax1.set_ylabel("Value")
    ax1.plot(x, y_obj, label="Objective f(x)", color="navy", linewidth=2, zorder=1)
    colors =['orange', 'green', 'purple', 'brown', 'cyan']

    if num_constraints > 0:
        ax1.axhline(0, color='black', linewidth=1, linestyle='-', label="Constraint Boundary (0)", zorder=0)
        for i, g_vals in enumerate(constraints):
            color = colors[i % len(colors)]
            ax1.plot(x, g_vals, linestyle="--", alpha=0.7, linewidth=1.5, label=f"Constraint g{i+1}(x)", color=color, zorder=1)

        y_min_fill = min(np.min(y_obj), -5)
        y_max_fill = max(np.max(y_obj), 5)
        max_constr_val = np.max([np.max(c) for c in constraints if np.max(c) < 1000]) if len(constraints) > 0 else 5
        y_max_fill = max(y_max_fill, max_constr_val)

        ax1.fill_between(x, -10000, 10000, where=infeasible_mask, color='red', alpha=0.1, label="Infeasible Region", zorder=0)
        ax1.set_ylim(bottom=np.min(y_obj) - 1, top=np.max(y_obj) + 1)
        current_ylim = ax1.get_ylim()
        ax1.set_ylim(min(current_ylim[0], -2), max(current_ylim[1], 2))

    if has_trials:
        feasible_trials_mask = (trial_v == num_constraints + 1)
        if np.any(feasible_trials_mask):
            ax1.scatter(trial_x[feasible_trials_mask], trial_z[feasible_trials_mask], color='blue', s=20, marker='o', alpha=0.6, label='Feasible Trials', zorder=3)

        infeasible_trials_mask = (trial_v <= num_constraints)
        if np.any(infeasible_trials_mask):
            ax1.scatter(trial_x[infeasible_trials_mask], trial_z[infeasible_trials_mask], color='red', s=20, marker='x', alpha=0.6, label='Infeasible Trials', zorder=3)

    ax1.scatter([found_x],[found_y], color="gold", s=150, marker='*', zorder=10, edgecolors='black', label="Found Minimum")
    ax1.legend(loc="best", fontsize=9, framealpha=0.9)
    ax1.grid(True, linestyle=':', alpha=0.6)

    ax2.set_title("Partially Computable Function View (Feasible Regions Only)")
    ax2.set_xlabel("X")
    ax2.set_ylabel("Objective Value")
    y_feasible_only = y_obj.copy()
    y_feasible_only[infeasible_mask] = np.nan
    ax2.plot(x, y_feasible_only, label="Feasible Objective f(x)", color="navy", linewidth=2, zorder=2)

    if has_trials:
        feasible_trials_mask = (trial_v == num_constraints + 1)
        if np.any(feasible_trials_mask):
            ax2.scatter(trial_x[feasible_trials_mask], trial_z[feasible_trials_mask], color='blue', s=25, marker='o', alpha=0.7, label='Feasible Trials', zorder=3)

    ax2.scatter([found_x], [found_y], color="gold", s=150, marker='*', zorder=10, edgecolors='black', label="Found Minimum")
    ax2.text(found_x, found_y, f"\n  Min\n  ({found_x:.3f}, {found_y:.3f})", verticalalignment='top', color='black', fontweight='bold', zorder=10)

    y_vals_clean = y_obj[~infeasible_mask]
    if len(y_vals_clean) > 0:
        ax2.set_ylim(np.min(y_vals_clean) - 0.5, np.max(y_vals_clean) + 0.5)

    ax2.legend(loc="best", fontsize=9, framealpha=0.9)
    ax2.grid(True, linestyle=':', alpha=0.6)

    plt.tight_layout()

    if silent_mode:
        os.makedirs("plots_output", exist_ok=True)
        plt.savefig(f"plots_output/1d_constr_task_{task_index}.png", dpi=150)
    else:
        plt.show()

if __name__ == "__main__":
    main()