import numpy as np 

from hadrana.ensembles.helpers import EnsembleHelpers
from hadrana.loader import load_rwfs, load_c2pt
from hadrana.momenta import get_momentum_shell
from hadrana.statistics import generate_jackknife_resamples


def compute_c2pt_jkn_fwd(ensemble: str, source_set: str, nsquare: int) -> np.ndarray:
    # compute averaged c2pt over measurements per nsquare value 

    c2pt_fwd: np.ndarray = load_c2pt(ensemble, source_set, "fwd")
    rwfs: np.ndarray     = load_rwfs(ensemble) 
    mom_shell_list, mom_shell_indices = get_momentum_shell(nsquare)

    print(mom_shell_list)
    print(mom_shell_indices)
    print(len(mom_shell_list))

    print()
    c2pt_shape = c2pt_fwd.shape
    c2pt_jkn_fwd = np.zeros(c2pt_shape, dtype=float)
    c2pt_jkn_fwd = generate_jackknife_resamples(c2pt_fwd[:, :, mom_shell_indices, :], rwfs, bin_size=1)

    print("SIZE BEFORE AVERAGING")
    print(f"{c2pt_jkn_fwd.nbytes / 1024 ** 2:.2f} MB")

    # average over source ids and equivalent momenta
    c2pt_jkn_fwd_avg = np.mean(c2pt_jkn_fwd, axis=(1,2))

    print("SIZE AFTER AVERAGING")
    print(f"{c2pt_jkn_fwd_avg.nbytes / 1024 ** 2:.2f} MB")

    return c2pt_jkn_fwd_avg

def compute_c2pt_jkn_fwd_bwd_avg(ensemble: str, nsquare: int) -> np.ndarray:
    ens = EnsembleHelpers(ensemble)

    source_sets: list[str] = ens.get_c2pt_source_sets()

    # load reweighting factors
    rwfs: np.ndarray = load_rwfs(ensemble) 

    # construct list of momentum configurations corresponding to nsquare
    mom_shell_list, mom_shell_indices = get_momentum_shell(nsquare)

    for source_set in source_sets:
        c2pt_fwd: np.ndarray = load_c2pt(ensemble, source_set, "fwd")
        c2pt_bwd: np.ndarray = load_c2pt(ensemble, source_set, "bwd")

        print(c2pt_fwd.shape)

    print(source_sets)


    pass

def compute_c2pt_ratio_jkn(ensemble: str, nsquare: int):
    pass

if __name__ == "__main__":
    ensemble = "D251"
    source_set = "source_set1"

    nsquare = 1

    c2pt_jkn_fwd_bwd_avg = compute_c2pt_jkn_fwd_bwd_avg(ensemble, nsquare)






