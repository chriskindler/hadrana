import numpy as np
import gvar as gv
import os
import pickle
import h5py as h5
from pathlib import Path

from hadrana.ensembles import *
from hadrana.loader import *

def get_trajectory(ensemble_id: str):
    ENSEMBLE_IDS_SYM = ["D250", "X251", "X250"]
    ENSEMBLE_IDS_TRM = ["D200", "D251", "N203"]
    ENSEMBLE_IDS_MSC = ["D201", "N201", "N204"]
    if ensemble_id in ENSEMBLE_IDS_SYM:
        trajectory = "SYM"
    elif ensemble_id in ENSEMBLE_IDS_TRM:
        trajectory = "TRM"
    elif ensemble_id in ENSEMBLE_IDS_MSC:
        trajectory = "MSC"
    else:
        raise ValueError(f"Unknown trajectory for ensemble {ensemble_id}.")
    return trajectory

def export_ratio_fit_results(ensemble_id: str, fit_result: dict, base_directory: str, overwrite: False):
    qmax = fit_result["QMAX"]
    fit_boundary_src = fit_result["FIT_RANGE"][0]
    fit_boundary_snk = fit_result["FIT_RANGE"][1]
    fit_model_id = fit_result["FIT_MODEL_ID"]

    filekey = f"{ensemble_id}-RATIO-FIT-RESULT-MODEL-{fit_model_id}-QMAX{qmax}-RANGE-SRC{fit_boundary_src}-SNK{fit_boundary_snk}"

    #TODO: Switch to h5 format for data pipeline
    # with h5.File(filename, "w") as f:
    #     group_p_initial = f.create_group("INITIAL_PARAMETERS")
    #     group_p_result  = f.create_group("RESULTING_PARAMETERS")
    #     group_matrix_element = f.create_group("MATRIX_ELEMENT")
    #     group_chi2 = f.create_group("CHI2")
    #     group_chi2dof = f.create_group("CHI2DOF")

    base_directory.mkdir(parents = True, exist_ok = True)
    filepath = base_directory / f"{filekey}.pkl"

    if overwrite:
        with open(filepath, "wb") as f:
            pickle.dump(fit_result, f)

def generate_form_factor_data(ensemble_ids: list[str], nsquares: list):
    # CONSTANTS
    hbarc = 0.1973269804

    # common lattice spacing for now
    a = gv.gvar("0.06379(37)")

    nbst = 100
    nmom = len(nsquares)
    nids = len(ensemble_ids)

    # matrix elements determined from ratio fits or sum ratio fits
    M_bst    = np.zeros((nbst, nids, nmom))

    # axial form factor data
    GA_bst = np.zeros((nbst, nids, nmom))

    # pseudo-bootstrap resamples
    a_bst = gv.sample(a, nbatch=nbst) 

    Mavg_bst = np.zeros((nbst, nids, nmom))
    Mspl_bst = np.zeros((nbst, nids, nmom))
    Qsquares = np.zeros((nbst, nids, nmom))

    for i, id in enumerate(ensemble_ids):
        # hadron masses
        MN = get_hadron_mass(id, HADRON_MASS_MAP, "nucleon")
        Mp = get_hadron_mass(id, HADRON_MASS_MAP, "pion")
        Mk = get_hadron_mass(id, HADRON_MASS_MAP, "kaon")

        # pseudo-bootstraps
        MN_bst = gv.sample(MN, nbatch=nbst)
        Mp_bst = gv.sample(Mp, nbatch=nbst)
        Mk_bst = gv.sample(Mk, nbatch=nbst)       

        N = get_spatial_extent(id)

        

        for j, n in enumerate(nsquares):
            qsquared = n * (2 * np.pi / N) ** 2
            EN_bst = np.sqrt(MN_bst ** 2 + qsquared)

            # kinematic factor E + M / sqrt{ 2E(E + M) }
            K_pbst = (EN_bst + MN_bst) / np.sqrt(2 * EN_bst * (EN_bst + MN_bst))

            if id == "D201":
                model_id = "A1"
                boundary_src = 3
                boundary_snk = 3
            else:
                boundary_src = 1
                boundary_snk = 1
                if n == 0:
                    model_id = "B3"
                else:
                    model_id = "D5"
            
            if id == "N201":
                if n == 0:
                    model_id = "B3"
                elif n in [1,2,4,5]:
                    model_id = "D5"
                    boundary_src = 1
                    boundary_snk = 1
                else:
                    model_id = "A1"
                    boundary_src = 4
                    boundary_snk = 4

            # load all bootstraps at once and renormalise them at bootstrap level
            M = load_ratio_fit_data(id, model_id, n, boundary_src, boundary_snk)["MATRIX_ELEMENT_BST"]

            # fill arrays
            M_bst[:, i, j]    = M 
            GA_bst[:, i, j] = M / K_pbst

            Mavg_bst[:, i, j] = (1 / 3) * (2 * Mk_bst ** 2 + Mp_bst ** 2)
            Mspl_bst[:, i, j] = Mk_bst ** 2 - Mp_bst ** 2

            """
            CONVERT TO GEV2
            """
            Qsquares[:, i, j] = qsquared * ((hbarc / a_bst) ** 2) - ((MN_bst - EN_bst) ** 2) * ((hbarc / a_bst) ** 2) 
    

    # STORE FORM FACTOR DATA PER ENSEMBLE
    form_factor_collection = {ensemble_id: {} for ensemble_id in ensemble_ids}

    for i, id in enumerate(ensemble_ids):
        GA_ens_cen = np.mean(GA_bst[:, i, :], axis=0)
        GA_ens_err = np.std(GA_bst[:, i, :], axis=0, ddof=1)
        GA_ens = gv.gvar(GA_ens_cen, GA_ens_err)

        Qsquares_ens_cen = np.mean(Qsquares[:, i, :], axis=0)

        form_factor_collection[id]["Q"] = Qsquares_ens_cen
        form_factor_collection[id]["GA"] = GA_ens

    return form_factor_collection


        
    

