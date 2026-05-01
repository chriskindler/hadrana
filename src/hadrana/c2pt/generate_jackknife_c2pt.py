import h5py
from hadrana.c2pt.compute_c2pt import compute_c2pt_jkn_fwd, compute_c2pt_jkn_fwd_bwd_avg, compute_c2pt_ratio_jkn

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

def export_c2pt_jkn(ensemble: str, nsquares: list[int]):
    export_path = f"/hdd/data/ensemble_data/{ensemble}/c2pt"
    export_name = f"{ensemble}_c2pt_jkn.h5"
    
    with h5py.File(f"{export_path}/{export_name}", "w") as f:
        f.attrs["ensemble"] = ensemble
        f.attrs["resampling"] = "jackknife"
        f.attrs["bin_size"] = 1
        f.attrs["reweighted"] = True
        f.attrs["exceptionals_excluded"] = True
        f.attrs["averaging"] = "source positions and equivalent momenta"
        for nsquare in nsquares:
            c2pt_fwd_jkn = compute_c2pt_jkn_fwd(ensemble, nsquare)
            f.create_dataset(f"c2pt/nsquare_{nsquare}/fwd", data = c2pt_fwd_jkn)

            c2pt_fwd_bwd_avg_jkn = compute_c2pt_jkn_fwd_bwd_avg(ensemble, nsquare)
            f.create_dataset(f"c2pt/nsquare_{nsquare}/fwd_bwd_avg", data = c2pt_fwd_bwd_avg_jkn)

            if nsquare != 0:
                c2pt_ratio_jkn = compute_c2pt_ratio_jkn(ensemble, nsquare)
                f.create_dataset(f"c2pt_ratio/nsquare_{nsquare}", data = c2pt_ratio_jkn)

if __name__ == "__main__":
    ensemble = "D251"
    nsquares = [0, 1, 2, 3, 4, 5, 6, 8]

    export_c2pt_jkn(ensemble, nsquares)    
