import h5py
import numpy as np
import pickle
from collections import defaultdict

from cora.correlators.twopoint import construct_c2pt_collection
from cora.ensembles.ensemble_helpers import *
from cora.reweighting import get_rwfs
from cora.statistics import generate_bootstrap_ensemble

"""
A650 run13_nog2 corresponds to source_set0
Source ids 0, 1, 2 correspond to sequential measurements using the LHPC trick
Source id 3 corresonds to a regular sequential measurement
"""

def dd():
    return defaultdict(dd)

def construct_isovector_c3pt(ensemble: str):
    ensemble = "A650"
    replicas = ['r000', 'r001']
    source_set = 'source_set0'
    source_ids  = [0, 1, 2, 3]
    t = 13

    Ncfg = get_Ncfg(ensemble)
    Nmsr = len(source_ids)
    Nmom = 93
    Ntau = t + 1

    # awful... split the function!
    C3A_isovector_per_replica  = dd()
    C3P_isovector_per_replica = dd()

    for replica in replicas:
        with h5py.File(f"/glurch/scratch/kic04594/cls_data/{ensemble}{replica}.h5") as f:
            ga = f["sequential3pt/current_channelAxialVector"]
            gp = f["sequential3pt/current_channelPseudoscalar"]

            gau = ga["current_flavorU"]; gad = ga["current_flavorD"]
            gpu = gp["current_flavorU"]; gpd = gp["current_flavorD"]

            # AXIAL 
            for i in ["x", "y", "z"]:
                di = f"current_direction{i}"
                gau_i = gau[di]; gad_i = gad[di]

                for j in ["X", "Y", "Z"]:
                    gau_data = gau_i[f"polarisation{j}/data"][source_set]
                    gad_data = gad_i[f"polarisation{j}/data"][source_set]

                    for s in source_ids:
                        mkey = str(s)
                        c3pt_U = gau_data[mkey]["fwd"][f"t{t}"][...]
                        c3pt_D = gad_data[mkey]["fwd"][f"t{t}"][...]

                        if np.isnan(c3pt_U).any() or np.isnan(c3pt_D).any():
                            print("NaNs in HDF5 raw data!",
                            "replica", replica, "cur", i, "pol", j, "s", s,
                            "U NaNs", np.isnan(c3pt_U).sum(),
                            "D NaNs", np.isnan(c3pt_D).sum())

                        val = c3pt_U - c3pt_D

                        C3A_isovector_per_replica[replica][i][j][t][s] = val
            
            # PSEUDOSCALAR
            for j in ["X", "Y", "Z"]:
                pol = f"polarisation{j}"
                if pol not in gpu:
                    continue
                gpu_data = gpu[f"{pol}/data"][source_set]
                gpd_data = gpd[f"{pol}/data"][source_set]

                for s in source_ids:
                    mkey = str(s)
                    c3pt_U = gpu_data[mkey]["fwd"][f"t{t}"][...]
                    c3pt_D = gpd_data[mkey]["fwd"][f"t{t}"][...]
                    val = c3pt_U - c3pt_D

                    C3P_isovector_per_replica[replica][j][t][s] = val

    # combine replicas per measurement
    C3A_isovector_per_meas = dd()
    for i in ['x', 'y', 'z']:
        for j in ['X', 'Y', 'Z']:
            for s in source_ids:
                C3A_isovector_per_meas[i][j][t][s] = np.concatenate(
                    [C3A_isovector_per_replica[r][i][j][t][s] for r in replicas], axis=0
                )

    C3P_isovector_per_meas = dd()
    for j in ['X', 'Y', 'Z']:
        for s in source_ids:
            C3P_isovector_per_meas[j][t][s] = np.concatenate(
                [C3P_isovector_per_replica[r][j][t][s] for r in replicas], axis=0
            )

    # build collections 
    C3A_isovector_collection = dd()
    for i in ['x', 'y', 'z']:
        for j in ['X', 'Y', 'Z']:
            arr = np.empty((Ncfg, Nmsr, Nmom, Ntau))
            for m, s in enumerate(source_ids):
                arr[:, m, :, :] = C3A_isovector_per_meas[i][j][t][s]
            C3A_isovector_collection[i][j][t] = arr

    C3P_isovector_collection = dd()
    for j in ['X', 'Y', 'Z']:
        arr = np.empty((Ncfg, Nmsr, Nmom, Ntau))
        for m, s in enumerate(source_ids):
            arr[:, m, :, :] = C3P_isovector_per_meas[j][t][s]
        C3P_isovector_collection[j][t]  = arr

    return (C3A_isovector_collection, C3P_isovector_collection)

def nan_report(name, x):
    x = np.asarray(x)
    print(f"{name}:")
    print("  shape:", x.shape)
    print("  NaNs:", np.isnan(x).sum())
    print("  Infs:", np.isinf(x).sum())
    print("  finite fraction:", np.isfinite(x).mean())
    print("  min/max:", np.nanmin(x), np.nanmax(x))
    print()

def construct_isovector_bootstrap_ratios(ensemble: str, cur: str, pol: str, Nbst: int, seed: int):
    # specific to A650
    t = 13
    source_ids  = [0, 1, 2, 3]

    # LOAD REWEIGHTING FACTORS 
    rwfs: np.ndarray           = get_rwfs(ensemble)
    nan_report("rwfs", rwfs)

    # LOAD TWO-POINT FUNCTIONS
    c2pt: np.ndarray  = construct_c2pt_collection(ensemble)
    # print(c2pt)
    print(c2pt.shape)

    # LOAD AXIAL AND PSEUDOSCALAR THREE-POINT FUNCTIONS
    c3pt_ax_coll, c3pt_ps_coll = construct_isovector_c3pt(ensemble)

    c3pt_ax = c3pt_ax_coll[cur][pol][t]

    c3pt_ps = c3pt_ps_coll[pol][t]

    Ncfg = c2pt.shape[0]

    assert rwfs.shape[0] == Ncfg
    assert c3pt_ax.shape[0] == Ncfg
    assert c3pt_ps.shape[0] == Ncfg

    # INIT RNG FOR BOOTSTRAPPING
    rng = np.random.default_rng(seed)
    bst_sample_idx = rng.integers(0, Ncfg, size=(Nbst, Ncfg), dtype=np.int64)
    
    # BOOTSTRAP C2PT
    c2pt_bst = generate_bootstrap_ensemble(c2pt, rwfs, Nbst, bst_sample_idx)

    # BOOTSTRAP C3PT
    c3pt_ax_bst = generate_bootstrap_ensemble(c3pt_ax, rwfs, Nbst, bst_sample_idx)
    c3pt_ps_bst = generate_bootstrap_ensemble(c3pt_ps, rwfs, Nbst, bst_sample_idx)

    # SANITY CHECKS
    print("\nTWO-POINT FUNCTIONS")
    print(c2pt_bst.shape)

    print("\nTHREE-POINT FUNCTIONS")
    print(c3pt_ax_bst.shape)
    print(c3pt_ps_bst.shape)

    # Nmsr = number of measurements = len([0, 1, 2, 3])
    Nmsr = len(source_ids) 
    Nmom = 93
    Ntau = t + 1

    """
    RATIO MAP 
    {
        "(13, 0)": 0, LHPC
        "(13, 1)": 1, LHPC
        "(13, 2)": 2, LHPC
        "(13, 3)": 3  SEQ
    }
    """
    
    # INIT RATIO NUMPY ARRAYS
    ratio_ax_bst = np.empty((Nbst, Nmsr, Nmom, Ntau))
    ratio_ps_bst = np.empty((Nbst, Nmsr, Nmom, Ntau))

    c2pt = c2pt_bst[:, :, 0, t][:, :, None, None]  # (Nbst, Nmsr, 1, 1)

    ratio_ax_bst = c3pt_ax_bst / c2pt
    ratio_ps_bst = c3pt_ps_bst / c2pt

    nan_report("ratio_ax_bst", ratio_ax_bst)
    nan_report("ratio_ps_bst", ratio_ps_bst)
    
    ratio_collection = dd()   

    ratio_collection['axial'] = ratio_ax_bst
    ratio_collection['pseudoscalar'] = ratio_ps_bst

    nan_report("c2pt (raw)", c2pt)
    nan_report("c3pt_ax (raw)", c3pt_ax)
    nan_report("c3pt_ps (raw)", c3pt_ps)

    return ratio_collection

def main():
    ensemble = "A650"
    Nbst = 100
    seed = 12345

    def to_dict(x):
        # for some reason defaultdicts cannot be pickled...
        # we need to convert to regular python dictionary
        if isinstance(x, defaultdict):
            x = {k: to_dict(v) for k, v in x.items()}
            return x
        if isinstance(x, dict):
            return {k: to_dict(v) for k, v in x.items()}
        return x

    ratios = construct_isovector_bootstrap_ratios(ensemble, 'z', 'Z', Nbst, seed)

    with open(f'/glurch/scratch/kic04594/data/misc/{ensemble}-LHPC-RATIOS-COMPARISON.pkl', "wb") as f:
        pickle.dump(to_dict(ratios), f)

if __name__  == "__main__":
    main()
