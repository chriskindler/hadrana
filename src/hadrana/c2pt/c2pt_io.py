import h5py
import numpy as np
import os
from pathlib import Path

from hadrana.ensembles.helpers import EnsembleHelpers
from hadrana.momenta import get_momentum_shell

def load_c2pt_raw(ensemble: str, nsquare: int, temporal_direction: str) -> np.ndarray:
    """
    Load c2pt correlators with shape (n_cfg, n_src, n_mom, n_time).
    only momentum configurations corresponding to momentum shell are read from disk
    we use 
        temporal_direction = "fwd" for c3pt/c2pt ratios
        average fwd/bwd for c2pt and c2pt ratio analysis
    """
    _, momentum_shell_indices = get_momentum_shell(nsquare)
    path = f"/hdd/data/ensemble_data/{ensemble}/c2pt/{ensemble}_cls2pt.h5"
    with h5py.File(path, "r") as f:
        c2pt = f[f"cls2pt/{temporal_direction}"][:, :, momentum_shell_indices, :]
    return c2pt

def load_c2pt_fwd_jkn(ensemble: str, source_set: str, nsquare: int):
    pass

def load_c2pt_fwd_bwd_avg_jkn(ensemble: str, nsquare: int):
    pass

def load_c2pt_ratio_jkn(ensemble: str, nsquare: int):
    pass



