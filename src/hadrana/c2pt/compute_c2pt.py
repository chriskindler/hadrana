import numpy as np 

from hadrana.c2pt.c2pt_io import load_c2pt_raw
from hadrana.ensembles.helpers import EnsembleHelpers
from hadrana.loader import load_rwfs
from hadrana.statistics import generate_jackknife_resamples

"""
D251_c2pt_jkn.h5/
    c2pt/
        nsquare_0/
            /fwd
            /fwd_bwd_avg
        nsquare_1/
            /fwd
            /fwd_bwd_avg
        nsquare_2/
            /fwd
            /fwd_bwd_avg
        nsquare_3/
            /fwd
            /fwd_bwd_avg
        nsquare_4/
            /fwd
            /fwd_bwd_avg
        nsquare_5/
            /fwd
            /fwd_bwd_avg
        nsquare_6/
            /fwd
            /fwd_bwd_avg
        nsquare_8/
            /fwd
            /fwd_bwd_avg

    c2pt_ratio/
        nsquare_1/
        nsquare_2/
        nsquare_3/
        nsquare_4/
        nsquare_5/
        nsquare_6/
        nsquare_8/
"""

def compute_c2pt_jkn_fwd(ensemble: str, nsquare: int, bin_size: int) -> np.ndarray:
    h = EnsembleHelpers(ensemble)
    exceptionals = h.get_exceptionals()

    rwfs: np.ndarray = load_rwfs(ensemble) 
    rwfs = np.delete(rwfs, exceptionals, axis=0)

    c2pt_fwd: np.ndarray = load_c2pt_raw(ensemble, nsquare, "fwd")
    c2pt_fwd = np.delete(c2pt_fwd, exceptionals, axis=0)

    c2pt_jkn_fwd = generate_jackknife_resamples(c2pt_fwd, rwfs, bin_size)

    # average over source ids and equivalent momenta
    c2pt_jkn_fwd_avg = np.mean(c2pt_jkn_fwd, axis=(1,2))
    return c2pt_jkn_fwd_avg

def compute_c2pt_jkn_fwd_bwd_avg(ensemble: str, nsquare: int, bin_size: int) -> np.ndarray:
    h = EnsembleHelpers(ensemble)
    exceptionals: np.ndarray = h.get_exceptionals()

    rwfs: np.ndarray = load_rwfs(ensemble) 
    rwfs = np.delete(rwfs, exceptionals, axis=0)

    c2pt_fwd: np.ndarray = load_c2pt_raw(ensemble, nsquare, "fwd")
    c2pt_bwd: np.ndarray = load_c2pt_raw(ensemble, nsquare, "bwd")

    c2pt_fwd = np.delete(c2pt_fwd, exceptionals, axis=0)
    c2pt_bwd = np.delete(c2pt_bwd, exceptionals, axis=0)

    # compute forward-backward average
    c2pt_fwd_bwd_avg = 0.5 * (c2pt_fwd + c2pt_bwd)  # (n_cfg, n_src, n_mom, n_time)

    # average over source ids and momenta
    c2pt = np.mean(c2pt_fwd_bwd_avg, axis=(1,2)) # (n_cfg, n_time)
    c2pt_jkn_fwd_bwd_avg = generate_jackknife_resamples(c2pt, rwfs, bin_size)
    return c2pt_jkn_fwd_bwd_avg

def compute_c2pt_ratio_jkn(ensemble: str, nsquare: int, bin_size: int):
    c2pt_jkn_q = compute_c2pt_jkn_fwd_bwd_avg(ensemble, nsquare, bin_size)
    c2pt_jkn_0 = compute_c2pt_jkn_fwd_bwd_avg(ensemble, 0, bin_size)
    return c2pt_jkn_q / c2pt_jkn_0

if __name__ == "__main__":
    ensemble = "D251"
    nsquare = 1
    bin_size = 2

    c2pt_jkn_fwd = compute_c2pt_jkn_fwd(ensemble, nsquare, bin_size)
    print(c2pt_jkn_fwd.shape)
    c2pt_jkn_fwd_bwd_avg = compute_c2pt_jkn_fwd_bwd_avg(ensemble, nsquare, bin_size)
    print(c2pt_jkn_fwd_bwd_avg.shape)
    c2pt_ratio_jkn = compute_c2pt_ratio_jkn(ensemble, nsquare, bin_size)
    print(c2pt_ratio_jkn.shape)



