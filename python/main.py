"""
Main entry point for the Satellite-UAV-UE Network Optimisation (Python).

Usage
-----
    python main.py                  # quick smoke-test
    python main.py --mode compare   # run clustering/caching comparison
    python main.py --mode converge  # run GA convergence study
    python main.py --mode all       # run both
"""
import argparse
import time

import numpy as np

import config as cfg
from model import Model
from initialize import create_model, initialize
from channel import path_sat_uav, path_uav_uav, path_uav_ue
from clustering import (game_theory_clustering, random_clustering,
                        nearest_clustering, random_cache)
from genetic_algorithm import run_ga
from compute_latency import compute_sum_latency
from deep_unfolding import deep_unfolding_cache


# ---------------------------------------------------------------------------
# Helper: build a model instance
# ---------------------------------------------------------------------------

def build_model(num_ue, num_ue_per_uav):
    """Construct and return an initialised Model."""
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
    return model, A, B


# ---------------------------------------------------------------------------
# Quick smoke-test
# ---------------------------------------------------------------------------

def smoke_test():
    """Run a single iteration to verify all modules work end-to-end."""
    print("=" * 60)
    print("Smoke Test – single scenario, single Monte Carlo realisation")
    print("=" * 60)
    num_ue = cfg.NUM_UAV * cfg.NUM_UE_PER_UAV
    num_ue_per_uav = cfg.NUM_UE_PER_UAV

    t0 = time.time()
    model, A_ini, B_ini = build_model(num_ue, num_ue_per_uav)
    print(f"  Model built in {time.time() - t0:.2f}s")
    print(model)

    kw = dict(
        area=cfg.AREA, Q=cfg.Q, noise_var=cfg.NOISE_VAR, BW=cfg.BW,
        num_uav_per_sat=cfg.NUM_UAV_PER_SAT, num_ue_per_uav=num_ue_per_uav,
        num_sat=cfg.NUM_SAT, num_uav=cfg.NUM_UAV, num_ue=num_ue,
        pmax_sat=cfg.PMAX_SAT, pmax_uav=cfg.PMAX_UAV,
    )

    # Game Theory Clustering
    t0 = time.time()
    A_gt = game_theory_clustering(
        model, A_ini, B_ini,
        cfg.NUM_SAT, cfg.NUM_UAV, num_ue,
        cfg.NUM_UAV_PER_SAT, num_ue_per_uav,
        cfg.N_U, cfg.PMAX_SAT, cfg.PMAX_UAV,
        cfg.AREA, cfg.Q, cfg.NOISE_VAR, cfg.BW,
    )
    print(f"  GT Clustering done in {time.time() - t0:.2f}s")

    # Compute latency (GTRP)
    B_rand = random_cache(cfg.NUM_UAV, cfg.NUM_M, cfg.NUM_F)
    t_gtrp, _ = compute_sum_latency(model, A_gt, B_rand, **kw)
    print(f"  GTRP total latency : {t_gtrp:.4f} s")

    # Compute latency (RCRP)
    A_rand = random_clustering(cfg.NUM_SAT, cfg.NUM_UAV, num_ue, cfg.N_U)
    B_rand = random_cache(cfg.NUM_UAV, cfg.NUM_M, cfg.NUM_F)
    t_rcrp, _ = compute_sum_latency(model, A_rand, B_rand, **kw)
    print(f"  RCRP total latency : {t_rcrp:.4f} s")

    # Compute latency (NCRP)
    A_near = nearest_clustering(model, cfg.NUM_SAT, cfg.NUM_UAV, num_ue, cfg.N_U)
    B_rand = random_cache(cfg.NUM_UAV, cfg.NUM_M, cfg.NUM_F)
    t_ncrp, _ = compute_sum_latency(model, A_near, B_rand, **kw)
    print(f"  NCRP total latency : {t_ncrp:.4f} s")

    # GA caching (sector 0)
    t0 = time.time()
    B_ga, bf = run_ga(
        cfg.GA_POP_SIZE, cfg.GA_PC, cfg.GA_PM, 0,
        model, A_gt, B_ini,
        cfg.AREA, cfg.Q, cfg.NOISE_VAR, cfg.BW,
        cfg.NUM_UAV_PER_SAT, num_ue_per_uav,
        cfg.NUM_M, cfg.NUM_F,
        cfg.NUM_SAT, cfg.NUM_UAV, num_ue,
        cfg.PMAX_SAT, cfg.PMAX_UAV,
        max_generation=5,  # quick test
    )
    t_ga, _ = compute_sum_latency(model, A_gt, B_ga, **kw)
    print(f"  GA (1 sector, 5 gen) done in {time.time() - t0:.2f}s; "
          f"latency={t_ga:.4f} s, best_fit={bf[-1]:.4f}")

    # Deep unfolding caching
    t0 = time.time()
    B_du = deep_unfolding_cache(
        model, A_gt, cfg.NUM_UAV, cfg.NUM_F, cfg.NUM_M,
        num_layers=cfg.DU_NUM_LAYERS, alpha=cfg.DU_ALPHA, beta=cfg.DU_BETA,
    )
    t_du, _ = compute_sum_latency(model, A_gt, B_du, **kw)
    print(f"  Deep unfolding cache done in {time.time() - t0:.2f}s; latency={t_du:.4f} s")

    print("\nSmoke test PASSED ✓")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Satellite-UAV-UE Network Optimisation"
    )
    parser.add_argument(
        "--mode",
        choices=["test", "compare", "converge", "all"],
        default="test",
        help=(
            "test    – quick smoke-test (default)\n"
            "compare – Monte Carlo comparison of strategies\n"
            "converge – GA convergence analysis\n"
            "all     – run compare + converge"
        ),
    )
    parser.add_argument("--monte", type=int, default=10,
                        help="Number of Monte Carlo realisations (compare mode).")
    args = parser.parse_args()

    if args.mode == "test":
        smoke_test()

    elif args.mode == "compare":
        from result_compare import run_comparison
        run_comparison(num_monte=args.monte)

    elif args.mode == "converge":
        from result_convergence import run_convergence
        run_convergence()

    elif args.mode == "all":
        from result_compare import run_comparison
        from result_convergence import run_convergence
        run_comparison(num_monte=args.monte)
        run_convergence()


if __name__ == "__main__":
    main()
