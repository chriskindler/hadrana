import h5py
from hadrana.c2pt.compute_c2pt import compute_c2pt_jkn_fwd, compute_c2pt_jkn_fwd_bwd_avg, compute_c2pt_ratio_jkn

"""
D251_c2pt_jkn.h5/
    c2pt/
        nsquare00/
            /fwd
            /fwd_bwd_avg
        nsquare01/
            /fwd
            /fwd_bwd_avg
        nsquare02/
            /fwd
            /fwd_bwd_avg
        nsquare03/
            /fwd
            /fwd_bwd_avg
        nsquare04/
            /fwd
            /fwd_bwd_avg
        nsquare05/
            /fwd
            /fwd_bwd_avg
        nsquare06/
            /fwd
            /fwd_bwd_avg
        nsquare08/
            /fwd
            /fwd_bwd_avg

    c2pt_ratio/
        nsquare01/
        nsquare02/
        nsquare03/
        nsquare04/
        nsquare05/
        nsquare06/
        nsquare08/
"""

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
            print(f"c2pt fwd")

            c2pt_fwd_jkn = compute_c2pt_jkn_fwd(ensemble, nsquare, bin_size)
            f.create_dataset(f"c2pt/nsquare{nsquare:02d}/fwd", data = c2pt_fwd_jkn)

            print(f"c2pt fwd bwd avg")
            c2pt_fwd_bwd_avg_jkn = compute_c2pt_jkn_fwd_bwd_avg(ensemble, nsquare, bin_size)
            f.create_dataset(f"c2pt/nsquare{nsquare:02d}/fwd_bwd_avg", data = c2pt_fwd_bwd_avg_jkn)

            print(f"c2pt ratio")
            if nsquare != 0:
                c2pt_ratio_jkn = compute_c2pt_ratio_jkn(ensemble, nsquare, bin_size)
                f.create_dataset(f"c2pt_ratio/nsquare{nsquare}", data = c2pt_ratio_jkn)
            print(f"nsquare = {nsquare} done")
            print()

if __name__ == "__main__":
    ensemble = "D251"
    bin_size = 2
    nsquares = [0, 1, 2, 3, 4, 5, 6, 8]

    export_c2pt_jkn(ensemble, nsquares, bin_size)    