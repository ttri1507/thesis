"""
Utility: compute the total system latency summed over all UEs.

Equivalent MATLAB file
----------------------
* func_ComputeSumLatency.m  → compute_sum_latency()
"""
import numpy as np

from power_allocation import power_allocation
from latency import func_latency


def compute_sum_latency(model, A, B,
                        area, Q, noise_var, BW,
                        num_uav_per_sat, num_ue_per_uav,
                        num_sat, num_uav, num_ue,
                        pmax_sat, pmax_uav):
    """
    Compute total latency and per-UE latency/case vectors.

    Parameters
    ----------
    (all parameters as described in power_allocation and func_latency)

    Returns
    -------
    t_total : float
        Sum of latencies across all UEs.
    T : np.ndarray float, shape (2, numUE)
        Row 0: latency per UE; Row 1: service case per UE.
    """
    T = np.zeros((2, num_ue))

    # Compute power allocation once for the current (A, B)
    P_sat, P_uav_uav, P_uav_ue, ASU, AUU = power_allocation(
        model, A, B,
        num_sat, num_uav, num_ue,
        num_uav_per_sat, num_ue_per_uav,
        pmax_sat, pmax_uav,
    )

    for k in range(num_ue):
        f = model.RP[k]
        tk, case_tk = func_latency(
            model, k, f, A, B,
            P_sat, P_uav_uav, P_uav_ue,
            area, Q, noise_var, BW,
            ASU, AUU,
        )
        T[0, k] = tk
        T[1, k] = case_tk

    t_total = T[0, :].sum()
    return t_total, T
