import h5py
import os
from pathlib import Path

from hadrana.c2pt.compute_c2pt import compute_c2pt_jkn_fwd, compute_c2pt_jkn_fwd_bwd_avg, compute_c2pt_ratio_jkn 
from hadrana.c2pt.c2pt_io import *



if __name__ == "__main__":
    ensemble = "D251"
    source_set = "source_set1"

    nsquare_values = [0, 1, 2, 3, 4, 5, 6, 8]

    for nsquare in nsquare_values:
        c2pt_jkn_fwd = compute_c2pt_jkn_fwd(ensemble, source_set, nsquare)

        print()







