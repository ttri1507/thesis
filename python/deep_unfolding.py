"""
Deep unfolding style cache optimisation for UAV cache placement.
"""
import numpy as np


def deep_unfolding_cache(model, A, num_uav, num_f, num_m,
                         num_layers=8, alpha=0.6, beta=0.2):
    """
    Build cache placement B (num_uav x num_f) using an unrolled iterative update.

    The update unfolds a fixed number of layers that combine:
      1) local request demand per UAV
      2) neighbour message passing from UAV-UAV channel affinity

    Returns
    -------
    B : np.ndarray int, shape (num_uav, num_f)
        Binary cache matrix that respects per-UAV capacity num_m.
    """
    # Demand matrix: requests seen by each UAV for each file.
    demand = np.zeros((num_uav, num_f), dtype=float)
    for k in range(model.RP.size):
        f = int(model.RP[k])
        u = int(np.argmax(A[:, k]))
        demand[u, f] += 1.0

    # Row-normalise demand for stable iterative updates.
    row_sum = demand.sum(axis=1, keepdims=True)
    row_sum = np.where(row_sum == 0.0, 1.0, row_sum)
    x = demand / row_sum

    # UAV neighbourhood affinity from channel magnitude.
    aff = np.abs(model.H_UU).astype(float)
    np.fill_diagonal(aff, 0.0)
    aff_sum = aff.sum(axis=1, keepdims=True)
    aff_sum = np.where(aff_sum == 0.0, 1.0, aff_sum)
    w = aff / aff_sum

    # Unrolled layers.
    for _ in range(max(1, int(num_layers))):
        msg = w @ x
        x = (1.0 - alpha) * x + alpha * demand / row_sum + beta * msg
        x = np.maximum(x, 0.0)

    # Project to binary cache decision with top-num_m per UAV.
    B = np.zeros((num_uav, num_f), dtype=int)
    k_keep = min(num_m, num_f)
    for u in range(num_uav):
        ids = np.argpartition(x[u], -k_keep)[-k_keep:]
        B[u, ids] = 1
    return B
