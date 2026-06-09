import numpy as np

from hadrana.c2pt.c2pt_io import load_c2pt_raw
from hadrana.loader import load_rwfs, load_c2pt_per_nsquare
from hadrana.statistics import identify_small_rwf_indices, identify_exceptionals

if __name__ == "__main__":
    ensemble        = "J501"
    momentum_shells = [0, 1, 2, 3, 4, 5, 6, 8]

    rwfs, rwfs_path = load_rwfs(ensemble)
    print(f"Loading {rwfs.shape[0]} reweighting factors from {rwfs_path}.")

    rwfs_flagged_indices = identify_small_rwf_indices(rwfs).flatten()
    nflag_rwfs = rwfs_flagged_indices.shape[0]


    print(f"Total {nflag_rwfs} flagged rwf(s), Config number(s) {rwfs_flagged_indices}")

    # for nsquare in momentum_shells: 
    #     print(f"nsquare = {nsquare}")
    #     c2pt_raw = load_c2pt_raw(ensemble, nsquare, "fwd")
    #
    #     if c2pt_raw.shape[0] != rwfs.shape[0]:
    #         raise RuntimeError(f"{rwfs.shape[0]} reweighting factors does not match {c2pt_raw.shape[0]} configurations.")
    #
    #     print(c2pt_raw.shape)
