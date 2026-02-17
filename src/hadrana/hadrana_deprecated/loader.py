import pickle

from hadrana.ensembles import *

def load_bootstrap_data(
        ensemble_id: str,
        renormalised: bool = False,
    ):

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
        raise ValueError(f"Unknown trajectory for ensemble id {ensemble_id}.")

    if renormalised:
        with open(f"/home/ck/phd/data/bootstrap_data/{trajectory}/{ensemble_id}/{ensemble_id}-BOOTSTRAP-DATA-REN-COLLECTION.pkl", "rb") as f:
            print(f"Loading bootstrap data collection: {ensemble_id}-BOOTSTRAP-DATA-REN-COLLECTION.pkl")
            ensemble_dict = pickle.load(f)
    else:
        with open(f"/home/ck/phd/data/bootstrap_data/{trajectory}/{ensemble_id}/{ensemble_id}-BOOTSTRAP-DATA-COLLECTION.pkl", "rb") as f:
            print(f"Loading bootstrap data collection: {ensemble_id}-BOOTSTRAP-DATA-COLLECTION.pkl")
            ensemble_dict = pickle.load(f)
    return ensemble_dict

def load_ratio_fit_data(
        ensemble_id: str,
        model_id: str,
        qmax: int,
        boundary_src: int,
        boundary_snk: int
    ):
    trajectory = get_trajectory(ensemble_id)
    filekey = f"{ensemble_id}-RATIO-FIT-RESULT-MODEL-{model_id}-QMAX{qmax}-RANGE-SRC{boundary_src}-SNK{boundary_snk}"

    with open(f"/home/ck/phd/data/fit_data/ratios/{trajectory}/{ensemble_id}/QMAX{qmax}/{filekey}.pkl", "rb") as f:
        #print(f"Loading fit data collection from: {filekey}.pkl")
        ratio_fit_collection = pickle.load(f)

    return ratio_fit_collection