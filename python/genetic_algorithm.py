"""
Genetic Algorithm for UAV cache placement optimisation.

Equivalent MATLAB files
-----------------------
* FUNC_GA.m                    → run_ga()
* func_GA_InitializePopulation → _initialize_population()
* func_GA_Evaluate             → _evaluate()
* func_GA_SelectParent         → _select_parent()
* func_GA_Crossover            → _crossover()
* func_GA_Mutation             → _mutation()
"""
import numpy as np

from power_allocation import power_allocation
from latency import func_latency
from config import INFEASIBLE_PENALTY


# ---------------------------------------------------------------------------
# GA operator: initialise population
# ---------------------------------------------------------------------------

def _initialize_population(pop_size, num_uav_per_sat, num_m, num_f):
    """
    Create a random initial population.

    Returns
    -------
    population : np.ndarray int, shape (num_uav_per_sat, num_f, pop_size)
        Each "individual" is a binary matrix (UAV × package).
    """
    population = np.zeros((num_uav_per_sat, num_f, pop_size), dtype=int)
    for p in range(pop_size):
        for u in range(num_uav_per_sat):
            n_cache = np.random.randint(1, num_m + 1)
            ids = np.random.choice(num_f, n_cache, replace=False)
            population[u, ids, p] = 1
    return population


# ---------------------------------------------------------------------------
# GA operator: evaluate fitness
# ---------------------------------------------------------------------------

def _evaluate(population, pop_size, model, A, B_full,
               area, Q, noise_var, BW, s,
               num_uav_per_sat, num_ue_per_uav, num_m,
               num_sat, num_uav, num_ue, pmax_sat, pmax_uav):
    """
    Compute the total latency (fitness) for each individual in the population.

    Returns
    -------
    fitness   : np.ndarray float, (pop_size,)
    fit_min   : float  Minimum fitness value.
    best_ele  : int    Index of the best individual.
    """
    num_ue_s = num_uav_per_sat * num_ue_per_uav
    ue_start = num_ue_s * s
    ue_end   = num_ue_s * (s + 1)
    uav_start = num_uav_per_sat * s
    uav_end   = num_uav_per_sat * (s + 1)

    fitness = np.zeros(pop_size)

    for par in range(pop_size):
        B_ind = population[:, :, par]

        # Feasibility check: no UAV caches more than num_m packages
        if np.any(B_ind.sum(axis=1) > num_m):
            fitness[par] = INFEASIBLE_PENALTY
            continue

        B_full_temp = B_full.copy()
        B_full_temp[uav_start:uav_end, :] = B_ind

        P_sat_t, P_uu_t, P_ue_t, ASU, AUU = power_allocation(
            model, A, B_full_temp,
            num_sat, num_uav, num_ue,
            num_uav_per_sat, num_ue_per_uav,
            pmax_sat, pmax_uav,
        )

        total_k = 0.0
        for k in range(ue_start, ue_end):
            f = model.RP[k]
            tk, _ = func_latency(
                model, k, f, A, B_full_temp,
                P_sat_t, P_uu_t, P_ue_t,
                area, Q, noise_var, BW,
                ASU, AUU,
            )
            total_k += tk
        fitness[par] = total_k

    best_ele = int(np.argmin(fitness))
    return fitness, fitness[best_ele], best_ele


# ---------------------------------------------------------------------------
# GA operator: roulette-wheel selection with elitism
# ---------------------------------------------------------------------------

def _select_parent(population, pop_size, fitness, best_ele):
    pop_out = np.zeros_like(population)
    # Eliminate infeasible individuals
    id_die  = np.where(fitness > INFEASIBLE_PENALTY * 0.9)[0]
    id_live = np.setdiff1d(np.arange(pop_size), id_die)
    fit_live = fitness[id_live]

    # Roulette-wheel (minimisation: invert fitness)
    fit_new = fit_live.max() - fit_live
    denom   = fit_new.sum()
    if denom == 0:
        prob = np.ones(len(id_live)) / len(id_live)
    else:
        prob = fit_new / denom
    cdf = np.concatenate([[0.0], np.cumsum(prob)])

    for par in range(pop_size):
        if par == best_ele:
            pop_out[:, :, par] = population[:, :, best_ele]
        else:
            r = np.random.rand()
            idx = np.searchsorted(cdf, r, side='right') - 1
            idx = min(idx, len(id_live) - 1)
            pop_out[:, :, par] = population[:, :, id_live[idx]]

    return pop_out


# ---------------------------------------------------------------------------
# GA operator: single-point crossover
# ---------------------------------------------------------------------------

def _crossover(population, pop_size, num_uav_per_sat, num_m, num_f, Pc, best_ele):
    child = np.zeros_like(population)
    for par1 in range(pop_size):
        if par1 == best_ele:
            child[:, :, par1] = population[:, :, par1]
            continue
        # Choose a distinct mate
        par2 = par1
        while par2 == par1:
            par2 = np.random.randint(pop_size)
        if np.random.rand() < Pc:
            cp = np.random.randint(1, num_f)   # cut point
            child[:, :cp, par1]   = population[:, :cp, par1]
            child[:, cp:, par1]   = population[:, cp:, par2]
        else:
            child[:, :, par1] = population[:, :, par1]
    return child


# ---------------------------------------------------------------------------
# GA operator: bit-flip mutation
# ---------------------------------------------------------------------------

def _mutation(population, Pm, best_ele):
    pop_out = population.copy()
    num_uav_per_sat, num_f, pop_size = population.shape
    n_bits = num_uav_per_sat * num_f

    for par in range(pop_size):
        if par == best_ele:
            continue
        if np.random.rand() < Pm:
            n_flip = max(1, int(np.ceil(n_bits * Pm)))
            ids    = np.random.choice(n_bits, n_flip, replace=False)
            for idx in ids:
                ui = idx % num_uav_per_sat
                fi = idx // num_uav_per_sat
                pop_out[ui, fi, par] ^= 1   # bit-flip
    return pop_out


# ---------------------------------------------------------------------------
# Main GA orchestrator
# ---------------------------------------------------------------------------

def run_ga(pop_size, Pc, Pm, s, model, A, B,
           area, Q, noise_var, BW,
           num_uav_per_sat, num_ue_per_uav, num_m, num_f,
           num_sat, num_uav, num_ue, pmax_sat, pmax_uav,
           max_generation=20, tolerance=5e-3, max_tole_gen=5):
    """
    Run the Genetic Algorithm to optimise cache placement for SAT sector s.

    Parameters
    ----------
    pop_size  : int    Population size.
    Pc        : float  Crossover probability.
    Pm        : float  Mutation probability.
    s         : int    Sector index (0-based).
    (remaining: standard system params)
    max_generation : int   Hard limit on generations.
    tolerance      : float Convergence threshold (seconds).
    max_tole_gen   : int   Allowed consecutive non-improving generations.

    Returns
    -------
    B_full   : np.ndarray int, (numUAV, numF)  Updated cache matrix.
    best_fit : list[float]                     Fitness vs. generation.
    """
    B_full = B.copy()
    uav_start = num_uav_per_sat * s
    uav_end   = num_uav_per_sat * (s + 1)

    # Initialise population
    parent = _initialize_population(pop_size, num_uav_per_sat, num_m, num_f)

    # Initial evaluation
    fitness, fit_min, best_ele = _evaluate(
        parent, pop_size, model, A, B_full,
        area, Q, noise_var, BW, s,
        num_uav_per_sat, num_ue_per_uav, num_m,
        num_sat, num_uav, num_ue, pmax_sat, pmax_uav,
    )
    best_fit = [fit_min]

    terminal = False
    generation = 0
    tole_generation = 0

    while not terminal:
        generation += 1

        parent = _select_parent(parent, pop_size, fitness, best_ele)
        parent = _crossover(parent, pop_size, num_uav_per_sat, num_m, num_f, Pc, best_ele)
        parent = _mutation(parent, Pm, best_ele)

        fitness, fit_min, best_ele = _evaluate(
            parent, pop_size, model, A, B_full,
            area, Q, noise_var, BW, s,
            num_uav_per_sat, num_ue_per_uav, num_m,
            num_sat, num_uav, num_ue, pmax_sat, pmax_uav,
        )
        best_fit.append(fit_min)

        if generation == max_generation:
            terminal = True
        elif generation > 1:
            if abs(fit_min - best_fit[-2]) < tolerance:
                tole_generation += 1
                if tole_generation == max_tole_gen:
                    terminal = True
            else:
                tole_generation = 0

    # Write best individual back to B_full
    B_full[uav_start:uav_end, :] = parent[:, :, best_ele]
    return B_full, best_fit
