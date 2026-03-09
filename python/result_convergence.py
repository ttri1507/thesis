"""
GA convergence study – compare different (P, Pc, Pm) configurations.

Equivalent MATLAB file
----------------------
* RESULT_Convergence.m  → run_convergence()
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import config as cfg
from model import Model
from initialize import create_model, initialize
from channel import path_sat_uav, path_uav_uav, path_uav_ue
from clustering import game_theory_clustering
from genetic_algorithm import run_ga


def run_convergence(
    output_file="result_convergence.png",
    verbose=True,
    pop_sizes=None,
    pc_list=None,
    pm_list=None,
):
    """
    Build a single model instance, run GA with multiple parameter configurations,
    and plot the convergence curve.

    Parameters
    ----------
    output_file : str   Path to save the result plot.
    verbose     : bool  Print progress.
    pop_sizes   : list[int]   Population sizes to try (default: [20, 40]).
    pc_list     : list[float] Crossover probabilities (default: [0.7, 0.8, 0.9]).
    pm_list     : list[float] Mutation probabilities  (default: [0.3, 0.2, 0.1]).
                  Each (pc, pm) pair is tried for every population size.

    Returns
    -------
    best_fits : list[list[float]]  Fitness-vs-generation for each config.
    legends   : list[str]          Legend strings for each config.
    """
    if pop_sizes is None:
        pop_sizes = [20, 40]
    if pc_list is None:
        pc_list = [0.7, 0.8, 0.9]
    if pm_list is None:
        pm_list = [0.3, 0.2, 0.1]
    num_ue = cfg.NUM_UAV * cfg.NUM_UE_PER_UAV

    # ---- Build model --------------------------------------------------------
    if verbose:
        print("Building network model …")
    model = Model()
    SAT, UAV, UE = create_model(cfg.AREA, cfg.NUM_SAT, cfg.NUM_UAV,
                                num_ue, cfg.Z_SAT, cfg.Z_UAV)
    model.SAT = SAT
    model.UAV = UAV
    model.UE  = UE

    model.H_SU = path_sat_uav(model, cfg.N)
    model.H_UU = path_uav_uav(model)
    model.H_UK = path_uav_ue(model, cfg.FC)

    H_norm2 = np.sqrt(np.abs(np.sum(np.conj(model.H_SU) * model.H_SU, axis=2)))
    H_norm2 = np.where(H_norm2 == 0, 1e-10, H_norm2)
    model.PrecodingS = np.conj(model.H_SU) / H_norm2[:, :, np.newaxis]

    RP, A, B, _, _, _ = initialize(
        cfg.NUM_SAT, cfg.NUM_UAV, num_ue, cfg.N_U,
        cfg.NUM_F, cfg.NUM_M, cfg.PMAX_SAT, cfg.PMAX_UAV,
    )
    model.RP = RP

    # ---- Game Theory Clustering (preprocessing) ----------------------------
    if verbose:
        print("Running Game Theory Clustering …")
    A = game_theory_clustering(
        model, A, B,
        cfg.NUM_SAT, cfg.NUM_UAV, num_ue,
        cfg.NUM_UAV_PER_SAT, cfg.NUM_UE_PER_UAV,
        cfg.N_U, cfg.PMAX_SAT, cfg.PMAX_UAV,
        cfg.AREA, cfg.Q, cfg.NOISE_VAR, cfg.BW,
    )

    # ---- GA parameter sweep ------------------------------------------------
    s = 0   # consider sector 0 (first SAT)

    best_fits = []
    legends   = []
    num_iterations = cfg.GA_MAX_GEN

    for P in pop_sizes:
        for Pc, Pm in zip(pc_list, pm_list):
            if verbose:
                print(f"  GA: P={P}, Pc={Pc}, Pm={Pm} …")
            _, bf = run_ga(
                P, Pc, Pm, s, model, A, B,
                cfg.AREA, cfg.Q, cfg.NOISE_VAR, cfg.BW,
                cfg.NUM_UAV_PER_SAT, cfg.NUM_UE_PER_UAV,
                cfg.NUM_M, cfg.NUM_F,
                cfg.NUM_SAT, cfg.NUM_UAV, num_ue,
                cfg.PMAX_SAT, cfg.PMAX_UAV,
                max_generation=num_iterations,
            )
            # Pad / truncate to exactly num_iterations points
            if len(bf) < num_iterations:
                bf = bf + [bf[-1]] * (num_iterations - len(bf))
            best_fits.append(bf[:num_iterations])
            legends.append(f"P={P}, Pc={Pc}, Pm={Pm}")

    # ---- Plot ---------------------------------------------------------------
    X = np.arange(1, num_iterations + 1)
    styles = ['r--^', 'r--o', 'r-*', 'b-^', 'b-o', 'b-*']
    fig, ax = plt.subplots(figsize=(10, 6))
    for i, (bf, leg) in enumerate(zip(best_fits, legends)):
        ax.semilogy(X, bf, styles[i % len(styles)], markersize=4, linewidth=1, label=leg)

    ax.set_xlabel('Generation')
    ax.set_ylabel('Best Fitness (Total Latency, s)')
    ax.set_title('GA Convergence for Different Parameter Configurations')
    ax.legend(fontsize=8)
    ax.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close(fig)
    if verbose:
        print(f"Plot saved to {output_file}")

    return best_fits, legends


if __name__ == "__main__":
    best_fits, legends = run_convergence()
    print("\nFinal fitness values per config:")
    for leg, bf in zip(legends, best_fits):
        print(f"  {leg}: {bf[-1]:.4f} s")
