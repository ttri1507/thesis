"""
Model data structure for the Satellite-UAV-UE network.
Equivalent to MATLAB's MODEL struct.
"""
import numpy as np


class Model:
    """
    Container for all network model data.

    Attributes
    ----------
    SAT : np.ndarray, shape (3, numSAT)
        3D positions [x, y, z] of satellites.
    UAV : np.ndarray, shape (3, numUAV)
        3D positions [x, y, z] of UAVs.
    UE : np.ndarray, shape (3, numUE)
        3D positions [x, y, z] of user equipment.
    RP : np.ndarray, shape (numUE,)
        Requested package index (0-based) for each UE.
    H_SU : np.ndarray, shape (numSAT, numUAV, N)
        Channel coefficients from SATs to UAVs (N antennas per SAT).
    H_UU : np.ndarray, shape (numUAV, numUAV)
        Channel coefficients between UAVs.
    H_UK : np.ndarray, shape (numUAV, numUE)
        Channel coefficients from UAVs to UEs.
    PrecodingS : np.ndarray, shape (numSAT, numUAV, N)
        Precoding vectors for each SAT-UAV link.
    IdClusUE : np.ndarray, shape (numUE,)
        Cluster (SAT index, 0-based) assignment for each UE.
    NumUE_in_SAT : np.ndarray, shape (numSAT,)
        Number of UEs assigned to each SAT cluster.
    """

    def __init__(self):
        self.SAT: np.ndarray = None
        self.UAV: np.ndarray = None
        self.UE: np.ndarray = None
        self.RP: np.ndarray = None
        self.H_SU: np.ndarray = None
        self.H_UU: np.ndarray = None
        self.H_UK: np.ndarray = None
        self.PrecodingS: np.ndarray = None
        self.IdClusUE: np.ndarray = None
        self.NumUE_in_SAT: np.ndarray = None

    def __repr__(self):
        parts = []
        for attr in ["SAT", "UAV", "UE", "RP", "H_SU", "H_UU", "H_UK", "PrecodingS"]:
            val = getattr(self, attr)
            if val is not None:
                parts.append(f"  {attr}: shape={val.shape}, dtype={val.dtype}")
            else:
                parts.append(f"  {attr}: None")
        return "Model(\n" + "\n".join(parts) + "\n)"
