"""
Power allocation for all three link types.

Equivalent MATLAB file
----------------------
* FUNC_PA.m  → power_allocation()
"""
import numpy as np


def power_allocation(model, A, B,
                     num_sat, num_uav, num_ue,
                     num_uav_per_sat, num_ue_per_uav,
                     pmax_sat, pmax_uav):
    """
    Water-filled (channel-aware weighted) power allocation.

    Parameters
    ----------
    model            : Model
    A                : np.ndarray int,   (numUAV, numUE)
    B                : np.ndarray int,   (numUAV, numF)
    num_sat          : int
    num_uav          : int
    num_ue           : int
    num_uav_per_sat  : int
    num_ue_per_uav   : int
    pmax_sat         : float
    pmax_uav         : float

    Returns
    -------
    P_sat      : np.ndarray float, (numSAT, numUAV)
    P_uav_uav  : np.ndarray float, (numUAV, numUAV)
    P_uav_ue   : np.ndarray float, (numUAV, numUE)
    ASU        : np.ndarray float, (numSAT, numUAV)   SAT→UAV request count
    AUU        : np.ndarray float, (numUAV, numUAV)   UAV→UAV request count
    """
    ASU = np.zeros((num_sat, num_uav))
    AUU = np.zeros((num_uav, num_uav))

    for k in range(num_ue):
        f  = model.RP[k]
        s  = k // (num_uav_per_sat * num_ue_per_uav)           # sector (0-based)
        uk_arr = np.where(A[:, k] == 1)[0]
        if len(uk_arr) == 0:
            continue
        uk = uk_arr[0]
        Us = np.arange(num_uav // num_sat * s, num_uav // num_sat * (s + 1))

        # Case 3: no UAV in sector has the package → SAT must serve
        ASU[s, uk] += int(not np.any(B[Us, f] == 1))

        # Case 2: primary UAV missing, but secondary UAV has it
        if (B[uk, f] == 0) and np.any(B[Us, f] == 1):
            Up = Us[B[Us, f] == 1]
            h_norms = np.abs(model.H_UU[Up, uk])
            best = np.argmax(h_norms)
            AUU[Up[best], uk] += 1

    # ---- SAT → UAV power ---------------------------------------------------
    # Effective channel power after beamforming: |h^H precoding|^2
    H_sq_su = np.abs(np.sum(np.conj(model.H_SU) * model.PrecodingS, axis=2)) ** 2  # (numSAT, numUAV)
    weight_s = np.where(H_sq_su > 0, ASU / H_sq_su, 0.0)
    denom_s  = weight_s.sum(axis=1, keepdims=True)
    denom_s  = np.where(denom_s == 0, 1.0, denom_s)  # avoid div-by-zero
    P_sat    = pmax_sat * weight_s / denom_s
    P_sat[weight_s == 0] = 0.0

    # ---- UAV → UAV power ---------------------------------------------------
    H_sq_uu = np.abs(model.H_UU) ** 2
    weight_u = np.where(H_sq_uu > 0, AUU / H_sq_uu, 0.0)
    denom_u  = weight_u.sum(axis=1, keepdims=True)
    denom_u  = np.where(denom_u == 0, 1.0, denom_u)
    P_uav_uav = pmax_uav * 0.25 * weight_u / denom_u
    P_uav_uav[weight_u == 0] = 0.0

    # ---- UAV → UE power ----------------------------------------------------
    H_sq_uk  = np.abs(model.H_UK) ** 2
    weight_uk = np.where(H_sq_uk > 0, A / H_sq_uk, 0.0)
    denom_uk  = weight_uk.sum(axis=1, keepdims=True)
    denom_uk  = np.where(denom_uk == 0, 1.0, denom_uk)
    P_uav_ue  = pmax_uav * 0.75 * weight_uk / denom_uk
    P_uav_ue[weight_uk == 0] = 0.0

    return P_sat, P_uav_uav, P_uav_ue, ASU, AUU
