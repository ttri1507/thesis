# Python Implementation – Satellite-UAV-UE Network Optimisation

This directory contains a complete Python rewrite of the original MATLAB codebase.

## System Overview

The system models a hybrid satellite-terrestrial network with:
- **4 Satellites (SATs)** at altitude 780 km  
- **16 UAVs** (4 per satellite) at altitude 0.5 km  
- **80 UEs** (ground users, configurable)  
- **30 content packages**, with each UAV caching up to 5  

The goal is to **minimise total latency** for content delivery via three approaches:
1. **Game Theory (GT) Clustering** – users greedily pick the best UAV in their sector  
2. **Genetic Algorithm (GA) Caching** – binary cache placement matrix evolved via GA  
3. Baselines: Random / Nearest-Neighbour clustering + Random caching  

---

## File Structure

| File | Description |
|------|-------------|
| `config.py` | System parameters (power, bandwidth, geometry, …) |
| `model.py` | `Model` data class (positions, channel matrices, etc.) |
| `initialize.py` | `create_model()`, `initialize()` – node placement & init |
| `channel.py` | Channel models: SAT→UAV (Shadowed-Rician), UAV↔UAV (Rayleigh), UAV→UE (ATG) |
| `latency.py` | Latency for all three service cases (direct / relay / satellite) |
| `power_allocation.py` | Weighted power allocation for SAT, UAV↔UAV, and UAV→UE links |
| `compute_latency.py` | `compute_sum_latency()` – aggregate latency over all UEs |
| `clustering.py` | GT clustering, random, nearest, coalition game, random cache |
| `genetic_algorithm.py` | Full GA pipeline (init, evaluate, select, crossover, mutate) |
| `plot_utils.py` | 2-D network topology visualisation |
| `result_compare.py` | Monte Carlo comparison of clustering/caching strategies |
| `result_convergence.py` | GA convergence study (different P, Pc, Pm) |
| `main.py` | CLI entry point |
| `requirements.txt` | Python dependencies |

---

## Installation

```bash
pip install -r requirements.txt
```

---

## Usage

```bash
# Quick smoke-test (verifies all modules work end-to-end)
python main.py

# Monte Carlo strategy comparison (default 10 realisations)
python main.py --mode compare --monte 10

# GA convergence study
python main.py --mode converge

# Run both comparison and convergence
python main.py --mode all
```

Output plots are saved as PNG files in the working directory.

---

## Key Differences from MATLAB

| Aspect | MATLAB | Python |
|--------|--------|--------|
| Indexing | 1-based | 0-based |
| Array storage | Column-major | Row-major (NumPy default) |
| Structures | `struct` | `Model` class |
| Random gamma | `gamrnd(a, b)` | `scipy.stats.gamma.rvs(a, scale=b)` |
| Matrix ops | `*` for element-wise | `*` in NumPy (also element-wise) |
| `find(…)` | returns indices | `np.where(…)[0]` |
| Complex math | `1i` | `1j` |

---

## Algorithm Details

### Channel Models
- **SAT→UAV**: Shadowed-Rician fading (ω=0.0005, σ²=0.063, m=2)  
- **UAV↔UAV**: Rayleigh fading with free-space path loss  
- **UAV→UE**: Air-to-Ground (ATG) model with LoS/NLoS probability  

### Service Cases (Latency)
1. Package cached at primary UAV → direct service  
2. Package at secondary UAV → relay via UAV-to-UAV link  
3. Package not cached → request from satellite  

### Genetic Algorithm
- **Chromosome**: Binary (numUAV_per_SAT × numF) cache matrix  
- **Fitness**: Total sector latency (lower = better)  
- **Selection**: Roulette-wheel with elitism  
- **Crossover**: Single-point, probability Pc  
- **Mutation**: Bit-flip, probability Pm  
