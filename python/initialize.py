"""
Network model creation and initialisation.

Equivalent MATLAB files
-----------------------
* FUNC_Create_Model.m
* FUNC_INITIALIZE.m
"""
import numpy as np

from model import Model


# ---------------------------------------------------------------------------
# FUNC_Create_Model
# ---------------------------------------------------------------------------

def create_model(area, num_sat, num_uav, num_ue, z_sat, z_uav):
    """
    Create 3-D positions for satellites, UAVs, and user equipment.

    Parameters
    ----------
    area : array-like, shape (2,)
        [x_area, y_area] in metres.
    num_sat : int
        Number of satellites (must be 4).
    num_uav : int
        Total number of UAVs.
    num_ue : int
        Total number of UEs.
    z_sat : float
        Satellite altitude (m).
    z_uav : float
        UAV altitude (m).

    Returns
    -------
    SAT : np.ndarray, shape (3, num_sat)
    UAV : np.ndarray, shape (3, num_uav)
    UE  : np.ndarray, shape (3, num_ue)
    """
    x_area, y_area = area[0], area[1]
    num_uav_per_sat = num_uav // num_sat
    num_ue_per_sat = num_ue // num_sat

    # ---- Satellites (fixed, assumed num_sat == 4) --------------------------
    SAT = np.zeros((3, 4))
    SAT[0, :] = [x_area / 4, x_area * 3 / 4, x_area * 3 / 4, x_area / 4]
    SAT[1, :] = [y_area / 4, y_area / 4, y_area * 3 / 4, y_area * 3 / 4]
    SAT[2, :] = z_sat

    # ---- UAVs (random within each SAT quadrant) ----------------------------
    uav_quadrants = [
        # (x_offset, y_offset)
        (0,           0),
        (x_area / 2,  0),
        (x_area / 2,  y_area / 2),
        (0,           y_area / 2),
    ]
    UAV_parts = []
    for x_off, y_off in uav_quadrants:
        block = np.full((3, num_uav_per_sat), z_uav)
        block[0, :] = x_off + x_area / 2 * np.random.rand(num_uav_per_sat)
        block[1, :] = y_off + y_area / 2 * np.random.rand(num_uav_per_sat)
        UAV_parts.append(block)
    UAV = np.concatenate(UAV_parts, axis=1)

    # ---- UEs (random within each SAT quadrant) -----------------------------
    UE_parts = []
    for x_off, y_off in uav_quadrants:
        block = np.zeros((3, num_ue_per_sat))
        block[0, :] = x_off + x_area / 2 * np.random.rand(num_ue_per_sat)
        block[1, :] = y_off + y_area / 2 * np.random.rand(num_ue_per_sat)
        UE_parts.append(block)
    UE = np.concatenate(UE_parts, axis=1)

    return SAT, UAV, UE


# ---------------------------------------------------------------------------
# FUNC_INITIALIZE
# ---------------------------------------------------------------------------

def initialize(num_sat, num_uav, num_ue, n_u, num_f, num_m, pmax_sat, pmax_uav):
    """
    Initialise requested packages (RP), user–UAV assignment (A),
    cache placement (B), and power matrices.

    Parameters
    ----------
    num_sat : int
    num_uav : int
    num_ue  : int
    n_u     : int     Max UEs per UAV.
    num_f   : int     Number of content packages.
    num_m   : int     Cache capacity per UAV.
    pmax_sat : float  SAT max power (W).
    pmax_uav : float  UAV max power (W).

    Returns
    -------
    RP          : np.ndarray int, shape (num_ue,)
                  Requested package index (0-based) for each UE.
    A           : np.ndarray int, shape (num_uav, num_ue)
                  Binary user-UAV assignment matrix.
    B           : np.ndarray int, shape (num_uav, num_f)
                  Binary cache placement matrix.
    P_sat       : np.ndarray float, shape (num_sat, num_uav)
    P_uav_uav   : np.ndarray float, shape (num_uav, num_uav)
    P_uav_ue    : np.ndarray float, shape (num_uav, num_ue)
    """
    num_uav_per_sat = num_uav // num_sat
    num_ue_per_sat = num_ue // num_sat

    # Requested packages (0-based)
    RP = np.random.randint(0, num_f, size=num_ue)

    # Cache placement: all UAVs cache the first num_m packages
    B = np.zeros((num_uav, num_f), dtype=int)
    B[:, :num_m] = 1

    # User-UAV assignment (diagonal round-robin per sector)
    A = np.zeros((num_uav, num_ue), dtype=int)
    for s in range(num_sat):
        uav_start = num_uav_per_sat * s
        ue_start = num_ue_per_sat * s
        for k in range(num_ue_per_sat):
            u = k % num_uav_per_sat
            A[uav_start + u, ue_start + k] = 1

    # SAT → UAV power (50 % of pmax_sat split equally across UAVs per SAT)
    P_sat = np.zeros((num_sat, num_uav))
    p_each_sat = 0.5 * pmax_sat / num_uav_per_sat
    for s in range(num_sat):
        col_start = num_uav_per_sat * s
        P_sat[s, col_start: col_start + num_uav_per_sat] = p_each_sat

    # UAV → UE power (75 % of pmax_uav, scaled by load)
    P_uav_ue = A * (pmax_uav * 0.75 / n_u)

    # UAV → UAV power (25 % of pmax_uav, uniform within sector, zero diagonal)
    P_uav_uav = np.zeros((num_uav, num_uav))
    p_sub_uu = pmax_uav * 0.25 / (num_uav_per_sat - 1)
    for s in range(num_sat):
        r = slice(num_uav_per_sat * s, num_uav_per_sat * (s + 1))
        block = np.full((num_uav_per_sat, num_uav_per_sat), p_sub_uu)
        np.fill_diagonal(block, 0.0)
        P_uav_uav[r, r] = block

    return RP, A, B, P_sat, P_uav_uav, P_uav_ue
