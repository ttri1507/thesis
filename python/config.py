"""
System configuration parameters.
All physical / scenario constants used across the codebase.
"""
import numpy as np

# ------------------------------------------------------------------
# Network topology
# ------------------------------------------------------------------
NUM_SAT = 4           # Number of satellites
NUM_UAV_PER_SAT = 4   # UAVs per satellite
NUM_UAV = NUM_UAV_PER_SAT * NUM_SAT   # Total UAVs  (16)
NUM_UE_PER_UAV = 10   # Average UEs per UAV
NUM_UE = NUM_UAV * NUM_UE_PER_UAV     # Total UEs   (80)

# ------------------------------------------------------------------
# Content / caching
# ------------------------------------------------------------------
NUM_F = 30   # Number of distinct content packages
NUM_M = 5    # Cache capacity per UAV (max packages stored)
N_U = 20     # Max UEs that one UAV can serve simultaneously

# ------------------------------------------------------------------
# Physical layer
# ------------------------------------------------------------------
N = 100            # Number of SAT antenna feeds / radiation elements
V_LIGHT = 3e8      # Speed of light  (m/s)
Q = 1e3            # Package size  (bits)
FC = 2e9           # Carrier frequency  (Hz)  – 2 GHz
BW = 20e6          # Bandwidth per link  (Hz)  – 20 MHz

# ------------------------------------------------------------------
# Power budgets
# ------------------------------------------------------------------
PMAX_SAT = 10 ** (50 / 10 - 3)   # SAT max power  ≈ 100 W  (50 dBm → W)
PMAX_UAV = 5.0                    # UAV max power   = 5 W

# ------------------------------------------------------------------
# Noise
# ------------------------------------------------------------------
NOISE_VARIANCE_DBM = -174 + 10 * np.log10(BW)          # Johnson–Nyquist noise
NOISE_VAR = 10 ** ((NOISE_VARIANCE_DBM - 30) / 10)     # Convert dBm → W

# ------------------------------------------------------------------
# Geometry
# ------------------------------------------------------------------
AREA = np.array([10e3, 10e3])   # Coverage area  10 km × 10 km
Z_SAT = 780e3                   # SAT altitude   780 km
Z_UAV = 0.5e3                   # UAV altitude   0.5 km

# ------------------------------------------------------------------
# GA parameters (defaults)
# ------------------------------------------------------------------
GA_POP_SIZE = 20    # Population size
GA_PC = 0.8         # Crossover probability
GA_PM = 0.2         # Mutation probability
GA_MAX_GEN = 20     # Maximum generations
GA_TOLERANCE = 5e-3 # Convergence tolerance (seconds)
GA_MAX_TOLE_GEN = 5 # Stop after this many consecutive non-improving generations

# ------------------------------------------------------------------
# Numerical constants
# ------------------------------------------------------------------
# Penalty fitness assigned to infeasible GA individuals
INFEASIBLE_PENALTY = 1e6
# Fallback latency returned when the computed value is non-finite or zero
FALLBACK_LATENCY = 1.0
