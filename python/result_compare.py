"""
Compare clustering + caching strategies via Monte Carlo simulation.

Equivalent MATLAB file
----------------------
* RESULT_Compare_ClusteringCaching.m  → run_comparison()

Six strategies compared:
  1. GTGA  – Game Theory clustering  + GA cache         (commented out – slow)
  2. RCGA  – Random clustering       + GA cache         (commented out – slow)
  3. NCGA  – Nearest clustering      + GA cache         (commented out – slow)
  4. GTRP  – Game Theory clustering  + Random cache
  5. RCRP  – Random clustering       + Random cache
  6. NCRP  – Nearest clustering      + Random cache
  7. DUGT  – Game Theory clustering  + Deep unfolding cache
"""
import numpy as np
import time
import matplotlib
matplotlib.use("Agg")   # non-interactive backend (safe for headless runs)
import matplotlib.pyplot as plt

import config as cfg
from model import Model
from initialize import create_model, initialize
from channel import path_sat_uav, path_uav_uav, path_uav_ue
from clustering import (game_theory_clustering, random_clustering,
                        nearest_clustering, random_cache)
from compute_latency import compute_sum_latency
from deep_unfolding import deep_unfolding_cache


def run_comparison(
    num_monte=10,
    num_ue_per_uav_list=None,
    output_file="result_compare.png",
    runtime_output_file="result_compare_runtime.png",
    verbose=True,
):
    """
    Run the Monte Carlo comparison experiment.

    Parameters
    ----------
    num_monte         : int   Number of Monte Carlo realisations.
    num_ue_per_uav_list : list[int]  Values of numUEperUAV to sweep.
    output_file       : str   Path to save the result plot.
    verbose           : bool  Print progress.

    Returns
    -------
    T_total : np.ndarray, dtype=float, shape (4, len(num_ue_per_uav_list))
              Average total latency for GTRP, RCRP, NCRP, DUGT strategies.
    R_total : np.ndarray, dtype=float, shape (4, len(num_ue_per_uav_list))
              Average execution time for GTRP, RCRP, NCRP, DUGT strategies.
    NumUE   : np.ndarray int
              Corresponding total UE counts.
    """
    if num_ue_per_uav_list is None:
        num_ue_per_uav_list = [5, 6, 7, 8, 9, 10]

    num_scenarios = len(num_ue_per_uav_list)
    num_strategies = 4   # GTRP, RCRP, NCRP, DUGT
    T_total = np.zeros((num_strategies, num_scenarios))
    R_total = np.zeros((num_strategies, num_scenarios))
    NumUE   = np.array([cfg.NUM_UAV * n for n in num_ue_per_uav_list])

    for scen_idx, num_ue_per_uav in enumerate(num_ue_per_uav_list):
        num_ue = cfg.NUM_UAV * num_ue_per_uav
        if verbose:
            print(f"Scenario {scen_idx + 1}/{num_scenarios}: numUE={num_ue}")

        T_monte = np.zeros((num_strategies, num_monte))
        R_monte = np.zeros((num_strategies, num_monte))

        for mc in range(num_monte):
            # ---- Build model -------------------------------------------------
            model = Model()
            SAT, UAV, UE = create_model(cfg.AREA, cfg.NUM_SAT, cfg.NUM_UAV,
                                        num_ue, cfg.Z_SAT, cfg.Z_UAV)
            model.SAT = SAT
            model.UAV = UAV
            model.UE  = UE

            model.H_SU = path_sat_uav(model, cfg.N)
            model.H_UU = path_uav_uav(model)
            model.H_UK = path_uav_ue(model, cfg.FC)

            # Conjugate beamforming precoding
            H_norm2 = np.sqrt(np.abs(np.sum(np.conj(model.H_SU) * model.H_SU, axis=2)))
            H_norm2 = np.where(H_norm2 == 0, 1e-10, H_norm2)
            model.PrecodingS = np.conj(model.H_SU) / H_norm2[:, :, np.newaxis]

            RP, A_ini, B_ini, _, _, _ = initialize(
                cfg.NUM_SAT, cfg.NUM_UAV, num_ue, cfg.N_U,
                cfg.NUM_F, cfg.NUM_M, cfg.PMAX_SAT, cfg.PMAX_UAV,
            )
            model.RP = RP

            kw = dict(
                area=cfg.AREA, Q=cfg.Q, noise_var=cfg.NOISE_VAR, BW=cfg.BW,
                num_uav_per_sat=cfg.NUM_UAV_PER_SAT, num_ue_per_uav=num_ue_per_uav,
                num_sat=cfg.NUM_SAT, num_uav=cfg.NUM_UAV, num_ue=num_ue,
                pmax_sat=cfg.PMAX_SAT, pmax_uav=cfg.PMAX_UAV,
            )

            # ---- Strategy 4: GTRP -------------------------------------------
            t_start = time.time()
            A_gt = game_theory_clustering(
                model, A_ini, B_ini,
                cfg.NUM_SAT, cfg.NUM_UAV, num_ue,
                cfg.NUM_UAV_PER_SAT, num_ue_per_uav,
                cfg.N_U, cfg.PMAX_SAT, cfg.PMAX_UAV,
                cfg.AREA, cfg.Q, cfg.NOISE_VAR, cfg.BW,
            )
            B_rand = random_cache(cfg.NUM_UAV, cfg.NUM_M, cfg.NUM_F)
            t_gtrp, _ = compute_sum_latency(model, A_gt, B_rand, **kw)
            r_gtrp = time.time() - t_start

            # ---- Strategy 5: RCRP -------------------------------------------
            t_start = time.time()
            A_rand = random_clustering(cfg.NUM_SAT, cfg.NUM_UAV, num_ue, cfg.N_U)
            B_rand = random_cache(cfg.NUM_UAV, cfg.NUM_M, cfg.NUM_F)
            t_rcrp, _ = compute_sum_latency(model, A_rand, B_rand, **kw)
            r_rcrp = time.time() - t_start

            # ---- Strategy 6: NCRP -------------------------------------------
            t_start = time.time()
            A_near = nearest_clustering(model, cfg.NUM_SAT, cfg.NUM_UAV, num_ue, cfg.N_U)
            B_rand = random_cache(cfg.NUM_UAV, cfg.NUM_M, cfg.NUM_F)
            t_ncrp, _ = compute_sum_latency(model, A_near, B_rand, **kw)
            r_ncrp = time.time() - t_start

            # ---- Strategy 7: DUGT -------------------------------------------
            t_start = time.time()
            A_du = game_theory_clustering(
                model, A_ini, B_ini,
                cfg.NUM_SAT, cfg.NUM_UAV, num_ue,
                cfg.NUM_UAV_PER_SAT, num_ue_per_uav,
                cfg.N_U, cfg.PMAX_SAT, cfg.PMAX_UAV,
                cfg.AREA, cfg.Q, cfg.NOISE_VAR, cfg.BW,
            )
            B_du = deep_unfolding_cache(
                model, A_du,
                cfg.NUM_UAV, cfg.NUM_F, cfg.NUM_M,
                num_layers=cfg.DU_NUM_LAYERS,
                alpha=cfg.DU_ALPHA,
                beta=cfg.DU_BETA,
            )
            t_dugt, _ = compute_sum_latency(model, A_du, B_du, **kw)
            r_dugt = time.time() - t_start

            T_monte[0, mc] = t_gtrp
            T_monte[1, mc] = t_rcrp
            T_monte[2, mc] = t_ncrp
            T_monte[3, mc] = t_dugt
            R_monte[0, mc] = r_gtrp
            R_monte[1, mc] = r_rcrp
            R_monte[2, mc] = r_ncrp
            R_monte[3, mc] = r_dugt

            if verbose:
                print(
                    f"  MC {mc + 1}/{num_monte}: "
                    f"t_gtrp={t_gtrp:.4f}, t_rcrp={t_rcrp:.4f}, "
                    f"t_ncrp={t_ncrp:.4f}, t_dugt={t_dugt:.4f} | "
                    f"r_gtrp={r_gtrp:.3f}s, r_rcrp={r_rcrp:.3f}s, "
                    f"r_ncrp={r_ncrp:.3f}s, r_dugt={r_dugt:.3f}s"
                )

        T_total[:, scen_idx] = T_monte.mean(axis=1)
        R_total[:, scen_idx] = R_monte.mean(axis=1)

    # ---- Plot ---------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(9, 6))
    styles = ['b-^', 'b-o', 'b-*', 'r-s']
    labels = ['GTRP', 'RCRP', 'NCRP', 'DUGT']
    for i in range(num_strategies):
        ax.semilogy(NumUE, T_total[i, :], styles[i], markersize=4, linewidth=1, label=labels[i])

    ax.set_xlabel('Number of UEs')
    ax.set_ylabel('Average Total Latency (s)')
    ax.set_title('Clustering & Caching Strategy Comparison')
    ax.legend()
    ax.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(output_file, dpi=150)
    plt.close(fig)
    if verbose:
        print(f"Plot saved to {output_file}")
        print("Average runtime across scenarios:")
        for i, label in enumerate(labels):
            print(f"  {label}: {R_total[i, :].mean():.3f}s")

    # Runtime plot
    fig, ax = plt.subplots(figsize=(9, 6))
    for i in range(num_strategies):
        ax.plot(NumUE, R_total[i, :], styles[i], markersize=4, linewidth=1, label=labels[i])
    ax.set_xlabel('Number of UEs')
    ax.set_ylabel('Average Execution Time (s)')
    ax.set_title('Strategy Runtime Comparison')
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(runtime_output_file, dpi=150)
    plt.close(fig)
    if verbose:
        print(f"Runtime plot saved to {runtime_output_file}")

    return T_total, R_total, NumUE


if __name__ == "__main__":
    T_total, R_total, NumUE = run_comparison(num_monte=3, num_ue_per_uav_list=[5, 7, 10])
    print("\nT_total (rows=strategies, cols=scenarios):")
    print(T_total)
    print("\nR_total (rows=strategies, cols=scenarios):")
    print(R_total)
