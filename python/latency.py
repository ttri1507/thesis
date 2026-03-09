"""
Latency computation for three service scenarios.

Equivalent MATLAB files
-----------------------
* FUNCLatency.m        → func_latency()
* funcLatencyUK.m      → latency_uk()
* funcLatencyUpK.m     → latency_upk()
* funcLatencySK.m      → latency_sk()
"""
import numpy as np

from config import FALLBACK_LATENCY


V_LIGHT = 3e8   # speed of light (m/s)


# ---------------------------------------------------------------------------
# Sub-functions (individual service cases)
# ---------------------------------------------------------------------------

def latency_uk(d_uk, h_uk, p_uk, h_other, p_other, Q, noise_var, BW):
    """
    Case 1 – UE served directly by its primary UAV.

    Parameters
    ----------
    d_uk    : float    Distance UAV→UE (m).
    h_uk    : complex  Channel coefficient UAV→UE.
    p_uk    : float    Transmit power UAV→UE (W).
    h_other : np.ndarray complex  Interfering channel coefficients on same UAV.
    p_other : np.ndarray float    Corresponding interfering powers.
    Q       : float    Package size (bits).
    noise_var : float  Noise variance (W).
    BW      : float    Bandwidth (Hz).

    Returns
    -------
    t_uk : float   Latency (s).
    """
    sinr = (np.abs(h_uk) ** 2 * p_uk) / (noise_var + np.sum(np.abs(h_other) ** 2 * p_other))
    r_uk = BW / np.log(2) * np.log1p(sinr)
    return 2 * d_uk / V_LIGHT + Q / r_uk


def latency_upk(d_uk, h_uk, p_uk, h_other, p_other,
                d_upu, h_upu, p_upu, Q, noise_var, BW, n_request):
    """
    Case 2 – UE served via a relay (secondary) UAV.

    Parameters
    ----------
    d_uk, h_uk, p_uk : primary UAV→UE link.
    h_other, p_other : interference on primary UAV.
    d_upu   : float    Distance secondary UAV → primary UAV (m).
    h_upu   : complex  Channel coefficient secondary→primary UAV.
    p_upu   : float    Transmit power secondary→primary UAV (W).
    n_request : int    Number of simultaneous requests on secondary UAV link.
    """
    sinr_uk = (np.abs(h_uk) ** 2 * p_uk) / (noise_var + np.sum(np.abs(h_other) ** 2 * p_other))
    r_uk    = BW / np.log(2) * np.log1p(sinr_uk)

    r_upu   = BW / np.log(2) * np.log1p(np.abs(h_upu) ** 2 * p_upu / noise_var)
    return 2 * (d_uk + d_upu) / V_LIGHT + Q / r_uk + Q * n_request / r_upu


def latency_sk(d_uk, h_uk, p_uk, h_other, p_other,
               d_su, h_su, precoding_su, p_su, Q, noise_var, BW, n_request):
    """
    Case 3 – UE served via satellite (no UAV caches the content).

    Parameters
    ----------
    d_su        : float     Distance SAT→UAV (m).
    h_su        : np.ndarray complex, shape (N,)  SAT→UAV channel vector.
    precoding_su: np.ndarray complex, shape (N,)  Precoding vector.
    p_su        : float     SAT transmit power (W).
    n_request   : int       Simultaneous requests on SAT→UAV link.
    """
    sinr_uk = (np.abs(h_uk) ** 2 * p_uk) / (noise_var + np.sum(np.abs(h_other) ** 2 * p_other))
    r_uk    = BW / np.log(2) * np.log1p(sinr_uk)

    # Beamformed effective channel: |h^H * w|^2
    effective = np.conj(h_su.ravel()) @ precoding_su.ravel()
    r_su = BW / np.log(2) * np.log1p(np.abs(effective) ** 2 * p_su / noise_var)
    return 2 * (d_uk + d_su) / V_LIGHT + Q / r_uk + Q * n_request / r_su


# ---------------------------------------------------------------------------
# Main dispatcher
# ---------------------------------------------------------------------------

def func_latency(model, k, f, A, B, P_sat, P_uav_uav, P_uav_ue,
                 area, Q, noise_var, BW, ASU, AUU):
    """
    Compute the latency for UE k requesting package f.

    Parameters
    ----------
    model  : Model
    k      : int    UE index (0-based).
    f      : int    Package index (0-based).
    A      : np.ndarray int,   shape (numUAV, numUE)   user-UAV assignment.
    B      : np.ndarray int,   shape (numUAV, numF)    cache placement.
    P_sat  : np.ndarray float, shape (numSAT, numUAV)
    P_uav_uav : np.ndarray float, shape (numUAV, numUAV)
    P_uav_ue  : np.ndarray float, shape (numUAV, numUE)
    area   : array-like [x_area, y_area].
    Q      : float  Package size (bits).
    noise_var : float
    BW     : float  Bandwidth (Hz).
    ASU    : np.ndarray float, shape (numSAT, numUAV)  SAT→UAV request counts.
    AUU    : np.ndarray float, shape (numUAV, numUAV)  UAV→UAV request counts.

    Returns
    -------
    tk      : float  Latency (s).
    case_tk : int    Service case (1, 2, or 3).
    """
    num_uav = model.UAV.shape[1]
    num_uav_per_sat = num_uav // model.SAT.shape[1]
    num_ue_per_uav = model.UE.shape[1] // num_uav

    # Sector index (0-based)
    s = k // (num_uav_per_sat * num_ue_per_uav)
    # UAV indices in this sector (0-based)
    us_slice = slice(num_uav // 4 * s, num_uav // 4 * (s + 1))
    Us = np.arange(num_uav // 4 * s, num_uav // 4 * (s + 1))

    # Primary UAV serving UE k
    u_arr = np.where(A[:, k] == 1)[0]
    if len(u_arr) == 0:
        return 1.0, 0
    u = u_arr[0]

    # Geometry & channel for UAV u – UE k link
    xu, yu, zu = model.UAV[:, u]
    xk, yk     = model.UE[0, k], model.UE[1, k]
    d_uk = np.sqrt((xu - xk) ** 2 + (yu - yk) ** 2 + zu ** 2)
    h_uk = model.H_UK[u, k]
    p_uk = P_uav_ue[u, k]

    # Interference from other UEs on the same primary UAV
    mask = np.ones(model.UE.shape[1], dtype=bool)
    mask[k] = False
    h_other = model.H_UK[u, mask]
    p_other = P_uav_ue[u, mask]

    # ---- Decide service case ------------------------------------------
    if B[u, f]:
        # Case 1: cached at primary UAV
        t_pri = latency_uk(d_uk, h_uk, p_uk, h_other, p_other, Q, noise_var, BW)
        tk = t_pri
        case_tk = 1

    elif np.any(B[Us, f] == 1):
        # Case 2: cached at a secondary UAV in the same sector
        Up = Us[B[Us, f] == 1]          # secondary UAVs that have the package
        h_norms = np.abs(model.H_UU[Up, u])
        best_local = np.argmax(h_norms)
        up = Up[best_local]

        xup, yup, zup = model.UAV[:, up]
        d_upu = np.sqrt((xup - xu) ** 2 + (yup - yu) ** 2 + (zup - zu) ** 2)
        h_upu = model.H_UU[up, u]
        p_upu = P_uav_uav[up, u]
        n_request = AUU[up, u]

        t_sec = latency_upk(d_uk, h_uk, p_uk, h_other, p_other,
                            d_upu, h_upu, p_upu, Q, noise_var, BW, n_request)
        tk = (1 - B[u, f]) * t_sec
        case_tk = 2

    else:
        # Case 3: request from satellite
        xs, ys, zs = model.SAT[:, s]
        d_su = np.sqrt((xs - xu) ** 2 + (ys - yu) ** 2 + (zs - zu) ** 2)
        h_su = model.H_SU[s, u, :]
        precoding_su = model.PrecodingS[s, u, :]
        p_su = P_sat[s, u]
        n_request = ASU[s, u]

        t_sat = latency_sk(d_uk, h_uk, p_uk, h_other, p_other,
                           d_su, h_su, precoding_su, p_su, Q, noise_var, BW, n_request)
        tk = (1 - B[u, f]) * (1 - int(np.any(B[Us, f]))) * t_sat
        case_tk = 3

    # Sanity check
    if not np.isfinite(tk) or tk == 0:
        tk = FALLBACK_LATENCY
    return tk, case_tk
