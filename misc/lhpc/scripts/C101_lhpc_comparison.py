import h5py
import numpy as np
import pickle
from collections import defaultdict

from cora.correlators.c2pt import construct_c2pt_collection
from cora.ensembles.ensemble_helpers import get_Ncfg
from cora.reweighting import get_rwfs
from cora.statistics import generate_bootstrap_ensemble

def dd():
    return defaultdict(dd)

def construct_isovector_c3pt(ensemble: str):
    if ensemble == "C101":
        replicas    = ['r014']
        source_sets = ['source_set3', 'source_set4']
        source_ids  = [0, 1, 2, 3]
        t = 14
    elif ensemble == "A650":
        replicas = ['r000', 'r001']
        source_sets = ['source_set0']
        source_ids  = [0, 1, 2, 3]
        t = 13

    # awful... split the function!
    C3A_reg_isovector_per_replica  = dd()
    C3P_reg_isovector_per_replica  = dd()
    C3A_lhpc_isovector_per_replica = dd()
    C3P_lhpc_isovector_per_replica = dd()

    for replica in replicas:
        with h5py.File(f"/glurch/scratch/kic04594/cls_data/{ensemble}{replica}.h5") as f:
            ga = f["sequential3pt/current_channelAxialVector"]
            gp = f["sequential3pt/current_channelPseudoscalar"]

            gau = ga["current_flavorU"]; gad = ga["current_flavorD"]
            gpu = gp["current_flavorU"]; gpd = gp["current_flavorD"]
            
            for i in ["x", "y", "z"]:
                di = f"current_direction{i}"
                gau_i = gau[di]; gad_i = gad[di]

                for j in ["X", "Y", "Z"]:
                    gau_data = gau_i[f"polarisation{j}/data"]
                    gad_data = gad_i[f"polarisation{j}/data"]

                    for source_set in source_sets:
                        gau_ss = gau_data[source_set]
                        gad_ss = gad_data[source_set]

                        for k in source_ids:
                            mkey = str(k)
                            c3pt_U = gau_ss[mkey]["fwd"]["t14"][...]
                            c3pt_D = gad_ss[mkey]["fwd"]["t14"][...]
                            val = c3pt_U - c3pt_D

                            if source_set == "source_set4":
                                C3A_reg_isovector_per_replica[replica][i][j][t][k] = val
                            else:
                                C3A_lhpc_isovector_per_replica[replica][i][j][t][k] = val

            for j in ["X", "Y", "Z"]:
                pol = f"polarisation{j}"
                if pol not in gpu:
                    continue
                gpu_data = gpu[f"{pol}/data"]
                gpd_data = gpd[f"{pol}/data"]

                for source_set in source_sets:
                    gpu_ss = gpu_data[source_set]
                    gpd_ss = gpd_data[source_set]

                    for k in source_ids:
                        mkey = str(k)
                        c3pt_U = gpu_ss[mkey]["fwd"]["t14"][...]
                        c3pt_D = gpd_ss[mkey]["fwd"]["t14"][...]
                        val = c3pt_U - c3pt_D

                        if source_set == "source_set4":
                            C3P_reg_isovector_per_replica[replica][j][t][k] = val
                        else:
                            C3P_lhpc_isovector_per_replica[replica][j][t][k] = val

    # get the shape for all the data (Ncfg, Nmsr, Nmom, Ntau)
    sample = C3A_reg_isovector_per_replica[replicas[0]]['z']['Z'][t][source_ids[0]]
    Ncfg, Nmom, Ntau = sample.shape
    Nmsr = len(source_ids)

    # combine replicas per measurement
    C3A_reg_isovector_per_meas = dd()
    C3A_lhpc_isovector_per_meas = dd()
    for i in ['x', 'y', 'z']:
        for j in ['X', 'Y', 'Z']:
            for k in source_ids:
                C3A_reg_isovector_per_meas[i][j][14][k] = np.concatenate(
                    [C3A_reg_isovector_per_replica[r][i][j][t][k] for r in replicas], axis=0
                )
                C3A_lhpc_isovector_per_meas[i][j][14][k] = np.concatenate(
                    [C3A_lhpc_isovector_per_replica[r][i][j][t][k] for r in replicas], axis=0
                )

    C3P_reg_isovector_per_meas = dd()
    C3P_lhpc_isovector_per_meas = dd()
    for j in ['X', 'Y', 'Z']:
        for k in source_ids:
            C3P_reg_isovector_per_meas[j][14][k] = np.concatenate(
                [C3P_reg_isovector_per_replica[r][j][t][k] for r in replicas], axis=0
            )
            C3P_lhpc_isovector_per_meas[j][14][k] = np.concatenate(
                [C3P_lhpc_isovector_per_replica[r][j][t][k] for r in replicas], axis=0
            )

    # build collections 
    C3A_reg_isovector_collection = dd()
    C3A_lhpc_isovector_collection = dd()
    for i in ['x', 'y', 'z']:
        for j in ['X', 'Y', 'Z']:
            arr_reg  = np.empty((Ncfg, Nmsr, Nmom, Ntau))
            arr_lhpc = np.empty((Ncfg, Nmsr, Nmom, Ntau))
            for m, k in enumerate(source_ids):
                arr_reg[:, m, :, :]  = C3A_reg_isovector_per_meas[i][j][14][k]
                arr_lhpc[:, m, :, :] = C3A_lhpc_isovector_per_meas[i][j][14][k]
            C3A_reg_isovector_collection[i][j][t]  = arr_reg
            C3A_lhpc_isovector_collection[i][j][t] = arr_lhpc

    C3P_reg_isovector_collection = dd()
    C3P_lhpc_isovector_collection = dd()
    for j in ['X', 'Y', 'Z']:
        arr_reg  = np.empty((Ncfg, Nmsr, Nmom, Ntau))
        arr_lhpc = np.empty((Ncfg, Nmsr, Nmom, Ntau))
        for m, k in enumerate(source_ids):
            arr_reg[:, m, :, :]  = C3P_reg_isovector_per_meas[j][t][k]
            arr_lhpc[:, m, :, :] = C3P_lhpc_isovector_per_meas[j][t][k]
        C3P_reg_isovector_collection[j][14]  = arr_reg
        C3P_lhpc_isovector_collection[j][14] = arr_lhpc

    return (C3A_reg_isovector_collection, C3P_reg_isovector_collection,
            C3A_lhpc_isovector_collection, C3P_lhpc_isovector_collection)

def construct_isovector_bootstrap_ratios(ensemble: str, cur: str, pol: str, Nbst: int, seed: int):
    """
        axial current direction cur in ['x', 'y', 'z'] 
        polarisation direction  pol in ['X', 'Y', 'Z']
    """

    if ensemble == "C101":
        Ncfg = 2000
    elif ensemble == "A650":
        Ncfg = 5062

    # LOAD REWEIGHTING FACTORS 
    rwfs: np.ndarray           = get_rwfs(ensemble)

    # LOAD TWO-POINT FUNCTIONS
    c2pt_reg: np.ndarray  = construct_c2pt_collection(ensemble)['source_set4']
    c2pt_lhpc: np.ndarray = construct_c2pt_collection(ensemble)['source_set3']

    # LOAD AXIAL AND PSEUDOSCALAR THREE-POINT FUNCTIONS
    c3pt_ax_reg_coll, c3pt_ps_reg_coll, c3pt_ax_lhpc_coll, c3pt_ps_lhpc_coll = construct_isovector_c3pt(ensemble)

    c3pt_ax_reg = c3pt_ax_reg_coll[cur][pol][14]
    c3pt_ax_lhpc = c3pt_ax_lhpc_coll[cur][pol][14]

    c3pt_ps_reg = c3pt_ps_reg_coll[pol][14]
    c3pt_ps_lhpc = c3pt_ps_lhpc_coll[pol][14]

    # sanity checks
    # print("C2PT")
    # print(c2pt_reg.shape)
    # print(c2pt_lphc.shape)
    #
    # print("AXIAL C3PT")
    # print(c3pt_ax_reg.shape)
    # print(c3pt_ax_lhpc.shape)
    #
    # print("PSEUDOSCLAR C3PT")
    # print(c3pt_ps_reg.shape)
    # print(c3pt_ps_lhpc.shape)

    # INIT RNG FOR BOOTSTRAPPING
    rng = np.random.default_rng(seed)
    bst_sample_idx = rng.integers(0, Ncfg, size=(Nbst, Ncfg), dtype=np.int64)
    
    # BOOTSTRAP C2PT
    c2pt_reg_bst = generate_bootstrap_ensemble(c2pt_reg, rwfs, Nbst, bst_sample_idx)
    c2pt_lhpc_bst = generate_bootstrap_ensemble(c2pt_lhpc, rwfs, Nbst, bst_sample_idx)

    # BOOTSTRAP C3PT
    c3pt_ax_reg_bst = generate_bootstrap_ensemble(c3pt_ax_reg, rwfs, Nbst, bst_sample_idx)
    c3pt_ax_lhpc_bst = generate_bootstrap_ensemble(c3pt_ax_lhpc, rwfs, Nbst, bst_sample_idx)
    c3pt_ps_reg_bst = generate_bootstrap_ensemble(c3pt_ps_reg, rwfs, Nbst, bst_sample_idx)
    c3pt_ps_lhpc_bst = generate_bootstrap_ensemble(c3pt_ps_lhpc, rwfs, Nbst, bst_sample_idx)

    # SANITY CHECKS
    print("\nTWO-POINT FUNCTIONS")
    print(c2pt_reg_bst.shape)
    print(c2pt_lhpc_bst.shape)

    print("\nTHREE-POINT FUNCTIONS")
    print(c3pt_ax_reg_bst.shape)
    print(c3pt_ax_lhpc_bst.shape)
    print(c3pt_ps_reg_bst.shape)
    print(c3pt_ps_lhpc_bst.shape)

    print("c2pt_reg_bst",  c2pt_reg_bst.shape)
    print("c3pt_ax_reg_bst", c3pt_ax_reg_bst.shape)
    print("c2pt_lhpc_bst", c2pt_lhpc_bst.shape)
    print("c3pt_ax_lhpc_bst", c3pt_ax_lhpc_bst.shape)

    # check that c3pt collections are distinct per (cur,pol)
    a1 = c3pt_ax_reg_coll['x']['X'][14]
    a2 = c3pt_ax_reg_coll['z']['Z'][14]
    print("same object?", a1 is a2) # should be False
    print("allclose?", np.allclose(a1, a2)) # should usually be False

    # Nmsr = number of measurements = len([0, 1, 2, 3])
    Nmsr = 4 
    Nmom = 93
    Ntau = 14 + 1
    """
    RATIO MAP 
    {
        "(14, 0)": 0,
        "(14, 1)": 1,
        "(14, 2)": 2,
        "(14, 3)": 3
    }
    """
    
    # INIT RATIO NUMPY ARRAYS
    ratio_ax_reg_bst  = np.empty((Nbst, Nmsr, Nmom, Ntau))
    ratio_ax_lhpc_bst = np.empty((Nbst, Nmsr, Nmom, Ntau))
    ratio_ps_reg_bst  = np.empty((Nbst, Nmsr, Nmom, Ntau))
    ratio_ps_lhpc_bst = np.empty((Nbst, Nmsr, Nmom, Ntau))

    c2pt_reg  = c2pt_reg_bst[:, :, 0, 14][:, :, None, None]  # (Nbst, Nmsr, 1, 1)
    c2pt_lhpc = c2pt_lhpc_bst[:, :, 0, 14][:, :, None, None] # (Nbst, Nmsr, 1, 1)

    ratio_ax_reg_bst  = c3pt_ax_reg_bst  / c2pt_reg
    ratio_ax_lhpc_bst = c3pt_ax_lhpc_bst / c2pt_lhpc
    ratio_ps_reg_bst  = c3pt_ps_reg_bst  / c2pt_reg
    ratio_ps_lhpc_bst = c3pt_ps_lhpc_bst / c2pt_lhpc
    
    ratio_collection = dd()   

    ratio_collection['axial']['reg'] = ratio_ax_reg_bst
    ratio_collection['axial']['lhpc'] = ratio_ax_lhpc_bst
    ratio_collection['pseudoscalar']['reg'] = ratio_ps_reg_bst
    ratio_collection['pseudoscalar']['lhpc'] = ratio_ps_lhpc_bst

    return ratio_collection

def main():
    ensemble = "C101"
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

    mom_id = 0
    src_id = 0
    ratio_ax_reg  = ratios['axial']['reg'][:,src_id,mom_id,:]
    ratio_ax_lhpc = ratios['axial']['lhpc'][:,src_id,mom_id,:]

    print(10 * "-", "REGULAR", 10 * "-")
    y_reg_est = np.mean(ratio_ax_reg, axis=0)
    y_reg_err = np.std(ratio_ax_reg, axis=0, ddof=1)
    print("CEN")
    print(y_reg_est)
    print("ERR")
    print(y_reg_err)

    print(10 * "-", "LHPC", 10 * "-")
    y_lhpc_est = np.mean(ratio_ax_lhpc, axis=0)
    y_lhpc_err = np.std(ratio_ax_lhpc, axis=0, ddof=1)
    print("CEN")
    print(y_lhpc_est)
    print("ERR")
    print(y_lhpc_err)

    # with open('/glurch/scratch/kic04594/data/misc/C101-LHPC-RATIOS-COMPARISON.pkl', "wb") as f:
    #     pickle.dump(to_dict(ratios), f)

if __name__  == "__main__":
    main()
