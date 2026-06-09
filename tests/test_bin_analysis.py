import numpy as np

from hadrana.c2pt.compute_c2pt import compute_c2pt_jkn_fwd_bwd_avg
from hadrana.loader import load_rwfs
from hadrana.statistics import generate_binsize_interval, maximum_binsize

if __name__ == "__main__":
    ensemble = "J501"
    nsquare  = 0
    binsizes =  [1, 2, 4, 8, 16, 32, 64]

    rwfs, rwfs_path = load_rwfs(ensemble)

    ncfg = rwfs.shape[0]
    smax = maximum_binsize(ncfg)
    binsizes.append(smax)

    print(smax)
    print(rwfs.shape)

    print(f"COMPUTE JACKKNIFE RESSAMPLES FOR BINSIZES \n{binsizes}")
    
    for s in binsizes:
        print(f"BINSIZE = {s}")
        jkn = compute_c2pt_jkn_fwd_bwd_avg(
            ensemble=ensemble,
            nsquare=nsquare, 
            bin_size=s
        )

        print(jkn.shape)

