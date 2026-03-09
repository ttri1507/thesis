"""
2-D network topology visualisation.

Equivalent MATLAB file
----------------------
* FUNC_Plot_model_2D.m  → plot_model_2d()
"""
import numpy as np
import matplotlib.pyplot as plt


def plot_model_2d(model, area, mode="normal"):
    """
    Plot the network topology.

    Parameters
    ----------
    model : Model
    area  : array-like  [x_area, y_area] in metres.
    mode  : str
        "normal"  – shows SATs, UAVs, and UEs.
        "cluster" – adds lines from UEs to their assigned SAT cluster.
    """
    x_area, y_area = area[0], area[1]
    SAT = model.SAT
    UAV = model.UAV
    UE  = model.UE

    fig, ax = plt.subplots(figsize=(8, 8))

    # Area border
    border_x = [0, x_area, x_area, 0, 0]
    border_y = [0, 0, y_area, y_area, 0]
    ax.plot(border_x, border_y, 'k-', linewidth=2)

    if mode == "normal":
        for m in range(SAT.shape[1]):
            ax.plot(SAT[0, m], SAT[1, m], 'r^',
                    markersize=12, linewidth=2, label='SAT' if m == 0 else '')
        for u in range(UAV.shape[1]):
            ax.plot(UAV[0, u], UAV[1, u], 'b^',
                    markersize=7, linewidth=1, label='UAV' if u == 0 else '')
        ax.plot(UE[0, :], UE[1, :], 'k*',
                markersize=5, linewidth=1, label='UE')
        ax.legend(['Border', 'SAT', 'UAV', 'UE'])

    elif mode == "cluster":
        if model.IdClusUE is None:
            raise ValueError("model.IdClusUE is not set – run clustering first.")
        colors = ['blue', 'green', 'orange', 'purple']
        ax.plot(UE[0, :], UE[1, :], 'b^', markersize=6, linewidth=1, label='UE')
        for u in range(UE.shape[1]):
            m = model.IdClusUE[u]
            ax.plot([UE[0, u], SAT[0, m]], [UE[1, u], SAT[1, m]],
                    color=colors[m % len(colors)], linewidth=0.8, alpha=0.5)
        for m in range(SAT.shape[1]):
            ax.plot(SAT[0, m], SAT[1, m], 'r^',
                    markersize=12, linewidth=2, label=f'SAT {m}')
        ax.legend()

    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.set_title(f'Network Topology ({mode} mode)')
    ax.grid(True)
    plt.tight_layout()
    return fig, ax
