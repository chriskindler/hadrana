import h5py
import numpy as np
from numpy import s_

from pathlib import Path

from hadrana.ensembles.helpers import EnsembleHelpers
from hadrana.momenta import get_momentum_shell

def load_rwfs(ensemble: str) -> np.ndarray:
    # if "rqcd" in ensemble:
    #     return np.ones(len(settings.get_total_config_list(ensemble)))
    ens = EnsembleHelpers(ensemble)
    replica_slices: dict = ens.get_replica_slices()

    n_cfg: int = ens.get_total_configuration_number()
    rwfs: np.ndarray = np.zeros(n_cfg)

    basic_path = Path("/hdd/data/sign_rwt_cls")
    for replica, slice_ in replica_slices.items():
        if ensemble == "S201" and replica == "r002":
            rwf_path = basic_path/f"rwt_ildg/cls_no_signs/{ensemble}{replica}.rwms.txt"
        elif ensemble == "U102" and replica == "r002":
            rwf_path = basic_path/f"rwt_ildg/cls_with_signs/{ensemble}{replica}.rwms.txt"
        else:
            # We prefer rwfs computed using deflation with signs
            rwf_path = basic_path/f"rwt_ildg/dfl_with_signs/{ensemble}{replica}.rwms.txt"

            # search in rwt_ildg.orig next
            if not rwf_path.exists():
                rwf_path = basic_path/f"rwt_ildg.orig/{ensemble}{replica}.rwms.txt"

            # search in prod...
            if not rwf_path.exists():
                rwf_path = basic_path/f"rwt_ildg/dfl_no_prod_with_signs/{ensemble}{replica}.rwms.txt"

            # if it doesn't exist we check for stochastic with signs 
            if not rwf_path.exists():
                rwf_path = basic_path/f"rwt_ildg/cls_with_signs/{ensemble}{replica}.rwms.txt"

            # if these do not exist we use dfl without signs
            if not rwf_path.exists():
                rwf_path = basic_path/f"rwt_ildg/dfl_no_signs/{ensemble}{replica}.rwms.txt"

            # search in prod...
            if not rwf_path.exists():
                rwf_path = basic_path/f"rwt_ildg/dfl_no_prod_no_signs/{ensemble}{replica}.rwms.txt"

            # and lastly we attempt stochastic without signs
            if not rwf_path.exists():
                rwf_path = basic_path/f"rwt_ildg/cls_no_signs/{ensemble}{replica}.rwms.txt"

        if not rwf_path.exists():
            raise RuntimeError("RWF not found")
      
        with open(rwf_path, 'r') as f:
            cfg_rwf_tmp = np.loadtxt(f)
            _rwf = cfg_rwf_tmp[:,1] 
            _rwf*= cfg_rwf_tmp[:,2] 

        rwfs[slice_] = _rwf
    
    return rwfs

def load_c2pt_per_nsquare(ensemble: str, nsquare_values: list[int], bin_size: int): 
    c2pt_per_nsquare = {}
    path = f"/hdd/data/ensemble_data/{ensemble}/c2pt/{ensemble}_c2pt_binsize{bin_size:02d}_jkn.h5"
    with h5py.File(path, "r") as file:
        for nsquare in nsquare_values:
            c2pt_per_nsquare[nsquare] = file[f"/c2pt/nsquare{nsquare:02d}/fwd_bwd_avg"][()]
    return c2pt_per_nsquare


def load_fits(run_dir, ensemble, nsquare, model_id, correlation_type, bin_size):
    signature = (f"{ensemble}-c2pt-nsquare{nsquare:02d}-binsize{bin_size:02d}"
                 f"-tmin*-tmax*-{model_id}-{correlation_type}.h5")
    fits = []
    for p in sorted(run_dir.glob(signature)):
        with h5py.File(p, "r") as f:
            d = {
                "tmin":    int(f["scalars/fit_range_min"][()]),
                "tmax":    int(f["scalars/fit_range_max"][()]),
                "chi2dof": float(f["scalars/chi2dof"][()]),
                "E0":      float(f["dicts/params_cen"].attrs["E0"]),
                "E0_err_jkn":  float(f["dicts/params_err"].attrs["E0"]),
                "E0_err_hess": float(f["dicts/params_err_hesse"].attrs["E0"]),
                "t_ext":         f["arrays/fit_range_ext"][()].flatten(),
                "y_fit_ext":     f["arrays/y_fit_cen_ext"][()].flatten(),
                "y_fit_err_ext": f["arrays/y_fit_err_ext"][()].flatten(),
                "Eeff_cen_ext":  f["arrays/Eeff_cen_ext"][()].flatten(),
                "Eeff_err_ext":  f["arrays/Eeff_err_ext"][()].flatten(),
            }
            d["E0_err"] = d["E0_err_jkn"]   # primary = jackknife
        fits.append(d)
    return fits