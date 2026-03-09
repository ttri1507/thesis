"""
Clustering algorithms.

Equivalent MATLAB files
-----------------------
* FUNC_GameTheory_Clustering.m  → game_theory_clustering()
* FUNC_RandomClustering.m       → random_clustering()
* FUNC_NearestClustering.m      → nearest_clustering()
* FUNC_CoalitionGame.m          → coalition_game()
* FUNC_RandomCache.m            → random_cache()
"""
import numpy as np

from power_allocation import power_allocation
from latency import func_latency

# Convergence thresholds for game-theory clustering
_GT_CONVERGE_FRAC = 0.99   # fraction of unchanged assignments to declare convergence
_GT_MAX_LOOPS = 100         # hard iteration limit


# ---------------------------------------------------------------------------
# Game Theory Clustering (user-centric best-response)
# ---------------------------------------------------------------------------

def game_theory_clustering(model, A_init, B,
                            num_sat, num_uav, num_ue,
                            num_uav_per_sat, num_ue_per_uav,
                            n_u, pmax_sat, pmax_uav,
                            area, Q, noise_var, BW):
    """
    Iterative best-response game: each UE is greedily reassigned to the UAV
    in its sector that minimises its own latency.

    Parameters
    ----------
    model          : Model
    A_init         : np.ndarray int, (numUAV, numUE)  Initial assignment.
    B              : np.ndarray int, (numUAV, numF)   Cache placement.
    n_u            : int   Max UEs per UAV.
    (remaining parameters: standard system params)

    Returns
    -------
    A_out : np.ndarray int, (numUAV, numUE)  Optimised assignment.
    """
    A_old = A_init.copy()
    num_loop = 0

    while True:
        num_loop += 1
        A_min = A_old.copy()

        for k in range(num_ue):
            f  = model.RP[k]
            s  = k // (num_uav_per_sat * num_ue_per_uav)
            Us = np.arange(num_uav // 4 * s, num_uav // 4 * (s + 1))

            ut_min = 1e6
            for id_u in range(num_uav_per_sat):
                u = Us[id_u]
                A_temp = A_min.copy()
                A_temp[:, k] = 0
                A_temp[u, k] = 1

                if np.all(A_temp.sum(axis=1) <= n_u):
                    P_sat_t, P_uu_t, P_ue_t, ASU, AUU = power_allocation(
                        model, A_temp, B,
                        num_sat, num_uav, num_ue,
                        num_uav_per_sat, num_ue_per_uav,
                        pmax_sat, pmax_uav,
                    )
                    ut, _ = func_latency(
                        model, k, f, A_temp, B,
                        P_sat_t, P_uu_t, P_ue_t,
                        area, Q, noise_var, BW,
                        ASU, AUU,
                    )
                else:
                    ut = 1e6

                if ut < ut_min:
                    ut_min = ut
                    A_min = A_temp.copy()

        # Convergence: > 99 % of assignments unchanged, or max iterations
        unchanged = np.sum(A_min == A_old)
        if (unchanged > _GT_CONVERGE_FRAC * num_uav * num_ue) or (num_loop == _GT_MAX_LOOPS):
            return A_min
        A_old = A_min.copy()


# ---------------------------------------------------------------------------
# Random Clustering
# ---------------------------------------------------------------------------

def random_clustering(num_sat, num_uav, num_ue, n_u):
    """
    Randomly assign UEs to UAVs within their sector, respecting n_u capacity.

    Returns
    -------
    A : np.ndarray int, (numUAV, numUE)
    """
    num_uav_per_sat = num_uav // num_sat
    num_ue_per_sat  = num_ue  // num_sat

    A = np.zeros((num_uav, num_ue), dtype=int)
    for s in range(num_sat):
        Us = np.arange(num_uav_per_sat * s, num_uav_per_sat * (s + 1))
        Ks = np.arange(num_ue_per_sat  * s, num_ue_per_sat  * (s + 1))
        for k_local in range(num_ue_per_sat):
            sums = A[Us, :].sum(axis=1)
            ok   = np.where(sums < n_u)[0]
            if len(ok) == 0:
                ok = np.arange(num_uav_per_sat)
            u = Us[np.random.choice(ok)]
            A[u, Ks[k_local]] = 1
    return A


# ---------------------------------------------------------------------------
# Nearest Clustering
# ---------------------------------------------------------------------------

def nearest_clustering(model, num_sat, num_uav, num_ue, n_u):
    """
    Assign each UE to the nearest non-overloaded UAV in its sector.

    Returns
    -------
    A : np.ndarray int, (numUAV, numUE)
    """
    num_uav_per_sat = num_uav // num_sat
    num_ue_per_sat  = num_ue  // num_sat

    # Euclidean distance (numUAV × numUE)
    diff = model.UAV[:, :, np.newaxis] - model.UE[:, np.newaxis, :]
    distance = np.sqrt(np.sum(diff ** 2, axis=0))

    A = np.zeros((num_uav, num_ue), dtype=int)
    for s in range(num_sat):
        Us = np.arange(num_uav_per_sat * s, num_uav_per_sat * (s + 1))
        Ks = np.arange(num_ue_per_sat  * s, num_ue_per_sat  * (s + 1))
        for k_local in range(num_ue_per_sat):
            k = Ks[k_local]
            sums  = A[Us, :].sum(axis=1)
            ok    = np.where(sums < n_u)[0]
            if len(ok) == 0:
                ok = np.arange(num_uav_per_sat)
            Us_ok = Us[ok]
            nearest = Us_ok[np.argmin(distance[Us_ok, k])]
            A[nearest, k] = 1
    return A


# ---------------------------------------------------------------------------
# Random Cache Placement
# ---------------------------------------------------------------------------

def random_cache(num_uav, num_m, num_f):
    """
    Randomly cache up to num_m packages at each UAV.

    Returns
    -------
    B : np.ndarray int, (numUAV, numF)
    """
    B = np.zeros((num_uav, num_f), dtype=int)
    for u in range(num_uav):
        n_cache = np.random.randint(1, num_m + 1)
        ids = np.random.choice(num_f, n_cache, replace=False)
        B[u, ids] = 1
    return B


# ---------------------------------------------------------------------------
# Coalition Game (SAT-level clustering)
# ---------------------------------------------------------------------------

def _utility(distance, chan, cluster, id_clus_ue, theta, t_total, h_total):
    """Utility value for a SAT cluster."""
    id_ue = np.where(id_clus_ue == cluster)[0]
    value = 0.0
    for u in id_ue:
        value += (theta * distance[cluster, u] / 3e8 / t_total
                  - (1 - theta) * chan[cluster, u] / h_total)
    return value


def coalition_game(model, B, t_thres, theta, area,
                   id_clus_ue_init, num_ue_in_sat_init):
    """
    Coalitional game at SAT-cluster level (join / swap operations).

    Parameters
    ----------
    model            : Model  (must have .Channel attribute for legacy code)
    B                : int   Max cluster size (beams).
    t_thres          : float Distance threshold in seconds (d_thres = t_thres * c).
    theta            : float Trade-off parameter ∈ [0, 1].
    area             : array-like
    id_clus_ue_init  : np.ndarray int, (numUE,)  Initial SAT assignment.
    num_ue_in_sat_init : np.ndarray int, (numSAT,)

    Returns
    -------
    model  (with updated .IdClusUE and .NumUE_in_SAT)
    """
    num_sat = model.SAT.shape[1]
    num_ue  = model.UE.shape[1]

    # Channel norm  (numSAT × numUE)
    chan = np.sqrt(np.abs(np.sum(np.conj(model.H_SU) * model.H_SU, axis=2)))
    h_total = chan.sum()

    # SAT-UE distance  (numSAT × numUE)
    diff = model.SAT[:, :, np.newaxis] - model.UE[:, np.newaxis, :]
    distance = np.sqrt(np.sum(diff ** 2, axis=0))
    t_total = distance.sum() / 3e8

    d_thres = t_thres * 3e8
    c_not   = distance > d_thres   # cannot-link matrix

    id_clus_ue     = id_clus_ue_init.copy()
    num_ue_in_sat  = num_ue_in_sat_init.copy()

    num_change_old = num_ue
    num_timechange = 0
    num_loop       = 0
    id_clus_ue_old = np.full(num_ue, -1)

    while not np.all(id_clus_ue_old == id_clus_ue):
        num_loop += 1
        id_clus_ue_old = id_clus_ue.copy()

        for u in range(num_ue):
            mp      = id_clus_ue[u]
            vmp_old = _utility(distance, chan, mp, id_clus_ue, theta, t_total, h_total)

            id_clus_ok    = id_clus_ue.copy()
            num_ue_sat_ok = num_ue_in_sat.copy()
            gap = 0.0

            # SATs that can be reached from UE u (no cannot-link), excluding mp
            m_ok = np.where(c_not[:, u] == 0)[0]
            m_ok = m_ok[m_ok != mp]

            for m in m_ok:
                id_clus_temp    = id_clus_ue.copy()
                num_ue_sat_temp = num_ue_in_sat.copy()
                vm_old = _utility(distance, chan, m, id_clus_ue, theta, t_total, h_total)

                if num_ue_sat_temp[m] == B:
                    # SWAP: exchange u with a member of cluster m
                    um = np.where(id_clus_ue == m)[0]
                    for up in um:
                        if c_not[mp, up]:
                            continue
                        id_clus_temp2 = id_clus_ue.copy()
                        id_clus_temp2[u]  = m
                        id_clus_temp2[up] = mp
                        vmp_new = _utility(distance, chan, mp, id_clus_temp2, theta, t_total, h_total)
                        vm_new  = _utility(distance, chan, m,  id_clus_temp2, theta, t_total, h_total)
                        gap_new = (vmp_old + vm_old) - (vmp_new + vm_new)
                        if gap_new > gap:
                            id_clus_ok = id_clus_temp2.copy()
                            gap = gap_new
                else:
                    # JOIN: move u from mp to m
                    id_clus_temp[u]  = m
                    num_ue_sat_temp[mp] -= 1
                    num_ue_sat_temp[m]  += 1
                    vmp_new = _utility(distance, chan, mp, id_clus_temp, theta, t_total, h_total)
                    vm_new  = _utility(distance, chan, m,  id_clus_temp, theta, t_total, h_total)
                    gap_new = (vmp_old + vm_old) - (vmp_new + vm_new)
                    if gap_new > gap:
                        id_clus_ok    = id_clus_temp.copy()
                        num_ue_sat_ok = num_ue_sat_temp.copy()
                        gap = gap_new

            id_clus_ue    = id_clus_ok.copy()
            num_ue_in_sat = num_ue_sat_ok.copy()

        num_change_new = int(np.sum(id_clus_ue_old != id_clus_ue))
        if (num_change_new - num_change_old) == 0:
            num_timechange += 1
        else:
            num_timechange = 0
        if num_timechange > 2 or num_loop == 100:
            break
        num_change_old = num_change_new

    model.IdClusUE      = id_clus_ue
    model.NumUE_in_SAT  = num_ue_in_sat
    return model
