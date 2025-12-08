import pickle

def load_bootstrap_data(
        id: str,
        renormalised: bool = False,
    ):
    if renormalised:
        with open(f"/home/ck/phd/data/bootstrap_data/{id}/{id}-BOOTSTRAP-DATA-REN-COLLECTION.pkl", "rb") as f:
            print(f"Loading bootstrap data collection: {id}-BOOTSTRAP-DATA-REN-COLLECTION.pkl")
            ensemble_dict = pickle.load(f)
    else:
        with open(f"/home/ck/phd/data/bootstrap_data/{id}/{id}-BOOTSTRAP-DATA-COLLECTION.pkl", "rb") as f:
            print(f"Loading bootstrap data collection: {id}-BOOTSTRAP-DATA-COLLECTION.pkl")
            ensemble_dict = pickle.load(f)
    return ensemble_dict

