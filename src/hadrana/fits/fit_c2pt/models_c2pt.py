# hadrana/fits/fit_c2pt/models.py
import numpy as np

# Convention: amplitudes A0, A1 are log-parametrized internally (Migrad sees logs);
# the gap dE1 is also log-parametrized (so E1 > E0 is guaranteed).

PARAMS_ONE_STATE = {"A0", "E0"}
PARAMS_TWO_STATE = {"A0", "E0", "A1", "dE1"}

def model_c2pt_one_state(x, A0, E0):
    return A0 * np.exp(-E0 * x)

def model_c2pt_two_state(x, A0, E0, A1, dE1):
    return A0 * np.exp(-E0 * x) + A1 * np.exp(-(E0 + dE1) * x)
MODEL_REGISTRY = {
    "one-state-exp": (model_c2pt_one_state, PARAMS_ONE_STATE),
    "two-state-exp": (model_c2pt_two_state, PARAMS_TWO_STATE),
}

def get_model(model_id: str):
    """Return (model_callable, set of log-parametrized parameters)."""
    if model_id not in MODEL_REGISTRY:
        raise ValueError(
            f"Unknown model_id={model_id!r}. Valid: {list(MODEL_REGISTRY)}"
        )
    return MODEL_REGISTRY[model_id]

def estimate_c2pt_starting_values(c2pt_jkn, t_start, t_final):
    c2pt_avg = np.mean(c2pt_jkn, axis=0)
    t_eff    = (t_start + t_final) // 2
    E0_start = float(np.log(c2pt_avg[t_eff] / c2pt_avg[t_eff + 1]))
    A0_start  = float(c2pt_avg[t_eff] * np.exp(E0_start * t_eff))
    return {
        "A0":  A0_start,
        "E0":  E0_start,
        "A1":  0.1 * A0_start,
        "dE1": 0.4,
    }


SCALAR_KEYS = (
    # identification
    "timestamp",
    "hash",
    
    "bin_size",
    "nsquare",

    "fit_id",
    "model_id",
    "resample_method",
    "correlation_type",

    # window
    "t_min",
    "t_max",

    # counts
    "n_data", # number of points included in the fit = len(np.arange(t_min, t_max + 1)) 
    "n_par",
    "n_dof",
    "n_res",
    "n_valid_res",

    # quality
    "chi2",
    "pvalue",
    "aic",
    "aicc",

    # Migrad flags
    "valid",
    "converged",
    "accurate", # return True if the covariance matrix is accuarate
    "fmin_edm",

    # Migrad inputs
    "tolerance",
    "strategy",
    "ncall",
)

ARRAY_KEYS = (
    "fit_range",       # (n_pts,) timeslices fit
    "fit_avg",         # (n_pts,) fit data central values
    "fit_err",         # (n_pts,) fit data errors
    "cov",             # (n_pts, n_pts)
    "cov_inv",         # (n_pts, n_pts)
    "residuals",       # (n_pts,) y_data - y_model at central values
    "valid_res",       # (n_res,) bool — convergence per resample
    "fit_range_ext",   # (n_ext,) extended grid for plotting
    "model_eval_ext",  # (n_res, n_dense) fit evaluated per resample on fine grid
)

DICT_KEYS = (
    "params_start",       # {name: float} - starting values
    "params_limits",      # {name: (lo, hi)} - bounds
    "params_res",         # {name: (n_res,) array} - resample values
    "params_cen",         # {name: float} - central values
    "params_err",         # {name: float} - jackknife/bootstrap error
    "params_err_hesse",   # {name: float} - Hesse error from migrad (cross-check with jackknife/bootstrap error)
)