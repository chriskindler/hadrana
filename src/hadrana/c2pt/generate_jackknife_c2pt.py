import h5py
import numpy as np
import time


from hadrana.c2pt.compute_c2pt import compute_c2pt_jkn_fwd, compute_c2pt_jkn_fwd_bwd_avg, compute_c2pt_ratio_jkn
from hadrana.statistics import maximum_binsize

def export_c2pt_jkn(ensemble: str, nsquares: list[int], bin_size: int):
    export_path = f"/hdd/data/ensemble_data/{ensemble}/c2pt"
    export_name = f"{ensemble}_c2pt_binsize{bin_size:02d}_jkn.h5"
    
    with h5py.File(f"{export_path}/{export_name}", "w") as f:
        f.attrs["ensemble"] = ensemble
        f.attrs["resampling"] = "jackknife"
        f.attrs["bin_size"] = bin_size
        f.attrs["reweighted"] = True
        f.attrs["exceptionals_excluded"] = True
        f.attrs["averaging"] = "source positions and equivalent momenta"
        print("Generate c2pt jackknife per nsquare")
        for nsquare in nsquares:
            print(f"nsquare = {nsquare}")
            c2pt_fwd_jkn = compute_c2pt_jkn_fwd(ensemble, nsquare, bin_size)
            print("Create c2pt_fwd dataset")
            f.create_dataset(f"c2pt/nsquare{nsquare:02d}/fwd", data = c2pt_fwd_jkn)

            c2pt_fwd_bwd_avg_jkn = compute_c2pt_jkn_fwd_bwd_avg(ensemble, nsquare, bin_size)
            print("Create c2pt_fwd_bwd_avg dataset")
            f.create_dataset(f"c2pt/nsquare{nsquare:02d}/fwd_bwd_avg", data = c2pt_fwd_bwd_avg_jkn)

            print(f"c2pt ratio")
            if nsquare != 0:
                c2pt_ratio_jkn = compute_c2pt_ratio_jkn(ensemble, nsquare, bin_size)
                f.create_dataset(f"c2pt_ratio/nsquare{nsquare}", data = c2pt_ratio_jkn)
            print(f"nsquare = {nsquare} done")
            print()

if __name__ == "__main__":
    ensemble = "A650"
    """
    smax(A650) = 50
    ncfg = 5062
    binsizes = [1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 30, 40, 50]

    smax(D251) = 20
    ncfg = 2012

    """
    momentum_shells = [0]

    from hadrana.loader import load_rwfs

    rwfs, rwfs_path = load_rwfs(ensemble)
    ncfg = rwfs.shape[0]

    binsizes = [1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 30, 40, maximum_binsize(ncfg)]

    print(f"ENSEMBLE: {ensemble}, CONFIGS = {ncfg}")
    print(f"BINSIZES = {binsizes}")
    print(f"MAXIMUM BINSIZE = {np.max(binsizes)}")

    timings = {}
    for s in binsizes:
        print(f"binsize = {s}")
        t0 = time.perf_counter()
        export_c2pt_jkn(ensemble, momentum_shells, s)
        dt = time.perf_counter() - t0
        timings[s] = dt
        print(f"binsize {s} took {dt:.2f}s")    

    print("\n TIMING SUMMARY ")
    for s, dt in timings.items():
        print(f"  binsize {s:>3}: {dt:7.2f}s")
    print(f"  total:     {sum(timings.values()):7.2f}s")
