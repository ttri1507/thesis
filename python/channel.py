"""
Channel / path-loss models.

Equivalent MATLAB files
-----------------------
* FUNC_Path_SAT_UAV.m   → path_sat_uav()
* FUNC_Path_UAV_UAV.m   → path_uav_uav()
* FUNC_Path_UAV_UE.m    → path_uav_ue()
"""
import numpy as np
from scipy.stats import gamma as gamma_dist


# ---------------------------------------------------------------------------
# SAT → UAV  (Shadowed-Rician fading)
# ---------------------------------------------------------------------------

def path_sat_uav(model, n_antennas):
    """
    Compute SAT-to-UAV channel coefficients using Shadowed-Rician fading.

    Parameters
    ----------
    model : Model
        Must have model.SAT (3, numSAT) and model.UAV (3, numUAV).
    n_antennas : int
        Number of radiation elements / antennas per SAT.

    Returns
    -------
    H : np.ndarray complex, shape (numSAT, numUAV, n_antennas)
    """
    num_sat = model.SAT.shape[1]
    num_uav = model.UAV.shape[1]

    # Euclidean distance matrix  (numSAT × numUAV)
    diff = model.SAT[:, :, np.newaxis] - model.UAV[:, np.newaxis, :]  # (3, numSAT, numUAV)
    distance = np.sqrt(np.sum(diff ** 2, axis=0))                      # (numSAT, numUAV)

    H = np.zeros((num_sat, num_uav, n_antennas), dtype=complex)
    for m in range(num_sat):
        for u in range(num_uav):
            H[m, u, :] = _channel_1sat_1uav(n_antennas, distance[m, u])
    return H


def _channel_1sat_1uav(n_antennas, dist):
    """Single SAT–UAV complex channel vector (Shadowed-Rician)."""
    pl_exponent = 2
    beta0 = 10 ** (-30 / 10)                          # reference path-loss gain
    large_scale = beta0 * dist ** (-pl_exponent)

    # Shadowed-Rician fading parameters
    om = 0.0005   # LoS component power  (ω)
    b  = 0.063    # scatter component    (σ²)
    m  = 2        # Nakagami-m parameter

    # Nakagami-m envelope
    naka_m = np.sqrt(gamma_dist.rvs(m, scale=1.0 / m))

    phi = 2 * np.pi * np.random.rand(n_antennas)
    p = np.sqrt(om) * np.cos(phi)
    q = np.sqrt(om) * np.sin(phi)

    sigma = np.sqrt(b)
    x = np.random.normal(naka_m * p, sigma)
    y = np.random.normal(naka_m * q, sigma)
    sr_fading = x + 1j * y

    return np.sqrt(large_scale) * sr_fading


# ---------------------------------------------------------------------------
# UAV → UAV  (Rayleigh fading)
# ---------------------------------------------------------------------------

def path_uav_uav(model):
    """
    Compute UAV-to-UAV channel coefficients (Rayleigh fading).

    Parameters
    ----------
    model : Model
        Must have model.UAV (3, numUAV).

    Returns
    -------
    H : np.ndarray complex, shape (numUAV, numUAV)
    """
    num_uav = model.UAV.shape[1]

    # Distance matrix  (numUAV × numUAV)
    diff = model.UAV[:, :, np.newaxis] - model.UAV[:, np.newaxis, :]
    distance = np.sqrt(np.sum(diff ** 2, axis=0))

    H = np.zeros((num_uav, num_uav), dtype=complex)
    for up in range(num_uav):
        for u in range(num_uav):
            H[up, u] = _channel_1uav_1uav(distance[up, u])
    return H


def _channel_1uav_1uav(dist):
    """Single UAV–UAV complex channel scalar (Rayleigh)."""
    pl_exponent = 2
    beta0 = 10 ** (-30 / 10)
    large_scale = beta0 * max(dist, 1e-6) ** (-pl_exponent)  # avoid div-by-zero
    small_scale = (np.random.randn() + 1j * np.random.randn()) / np.sqrt(2)
    return np.sqrt(large_scale) * small_scale


# ---------------------------------------------------------------------------
# UAV → UE  (Air-to-Ground path loss + Rayleigh fading)
# ---------------------------------------------------------------------------

def path_uav_ue(model, fc):
    """
    Compute UAV-to-UE channel coefficients (ATG path loss + Rayleigh).

    Parameters
    ----------
    model : Model
        Must have model.UAV (3, numUAV) and model.UE (3, numUE).
    fc : float
        Carrier frequency (Hz).

    Returns
    -------
    H : np.ndarray complex, shape (numUAV, numUE)
    """
    num_uav = model.UAV.shape[1]
    num_ue  = model.UE.shape[1]

    # 3-D Euclidean distance  (numUAV × numUE)
    diff = model.UAV[:, :, np.newaxis] - model.UE[:, np.newaxis, :]
    distance = np.sqrt(np.sum(diff ** 2, axis=0))

    H = np.zeros((num_uav, num_ue), dtype=complex)
    for u in range(num_uav):
        z_u = model.UAV[2, u]
        for k in range(num_ue):
            dist = distance[u, k]
            dist_xy = np.sqrt(max(dist ** 2 - z_u ** 2, 0.0))
            H[u, k] = _channel_1uav_1ue(dist, dist_xy, z_u, fc)
    return H


def _channel_1uav_1ue(dist, dist_xy, z_u, fc):
    """Single UAV–UE complex channel scalar (ATG + Rayleigh)."""
    pl_exponent = 2
    c = 3e8
    a_atg = 9.61
    b_atg = 0.16
    eta_los  = 1
    eta_nlos = 20

    pl_db = 10 * pl_exponent * np.log10(4 * np.pi * fc * max(dist, 1e-6) / c)

    if dist_xy < 1e-6:
        angle = np.pi / 2
    else:
        angle = np.arctan(z_u / dist_xy)

    pr_los  = 1.0 / (1.0 + a_atg * np.exp(-b_atg * (np.degrees(angle) - a_atg)))
    pr_nlos = 1.0 - pr_los

    path_loss_db = pl_db + eta_los * pr_los + eta_nlos * pr_nlos
    path_loss_lin = 10 ** (-path_loss_db / 10)

    small_scale = (np.random.randn() + 1j * np.random.randn()) / np.sqrt(2)
    return np.sqrt(path_loss_lin) * small_scale
