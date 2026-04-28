import h5py
import numpy as np

from hadrana.ensembles.helpers import EnsembleHelpers
from hadrana.momenta import get_momentum_list

def load_rwfs(ensemble: str) -> np.ndarray:
    """Load reweighting factors (RWTM2_EO * RWRAT) concatenated across replicas."""
    ens = EnsembleHelpers(ensemble)
    base = f"/hdd/data/ensemble_data/{ensemble}/rwfs"

    rwfs_per_replica = []
    for replica in ens.get_replicas():
        path = f"{base}/{ensemble}{replica}.rwms.txt"
        data = np.loadtxt(path)
        rwfs_per_replica.append(data[:, 1] * data[:, 2])

    return np.concatenate(rwfs_per_replica)

def load_c2pt(ensemble: str, source_set: str, temporal_direction: str) -> np.ndarray:
    """
    Load c2pt correlators with shape (n_cfg, n_src, n_mom, n_time).
    we use 
        temporal_direction = "fwd" for c3pt/c2pt ratios
        average fwd/bwd for c2pt and c2pt_ratio analysis
    """

    ens = EnsembleHelpers(ensemble)
    src_ids = ens.get_c2pt_measurement_ids(source_set)

    c2pt_shape = (
        ens.get_total_configuration_number(),
        len(src_ids),
        len(get_momentum_list()),
        ens.get_temporal_lattice_dimension(),
    )
    c2pt = np.zeros(c2pt_shape, dtype=float)
    path = f"/hdd/data/ensemble_data/{ensemble}/c2pt/{ensemble}_c2pt.h5"
    with h5py.File(path, "r") as f:
        for s, src_id in enumerate(src_ids):
            key = f"C2pt/polarisationNone/data/{source_set}/{src_id}/{temporal_direction}"
            c2pt[:, s, :, :] = f[key][()]

    return c2pt

def load_c3pt(ensemble: str, source_set: str):
    pass

if __name__ == "__main__":
    ensemble = "D251"

    source_set = "source_set1"
    temporal_direction = "fwd"

    c2pt_data = load_c2pt(ensemble)
