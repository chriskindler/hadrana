# hadrana/fits/results.py
from __future__ import annotations
import numpy as np
from scipy.stats import chi2 as chi2_dist


# The fit-result dict has three categories of values, each goes to a different
# place during serialization:
#   - SCALAR_KEYS  → manifest (parquet) and HDF5 attributes
#   - ARRAY_KEYS   → HDF5 datasets (too bulky for parquet)
#   - DICT_KEYS    → flattened into manifest columns; saved as nested groups in HDF5

SCALAR_KEYS = (
    # identification
    "fit_id",
    "model_id",
    "resample_method",
    "correlated",
    "bin_size",

    # window
    "t_min",
    "t_max",
    "n_pts", # number of points included in the fit

    # counts
    "n_par",
    "n_dof",
    "n_data",
    "n_res",
    "n_valid_res",

    # quality
    "chi2",
    "chi2_red",
    "pvalue",
    "aic",
    "bic",

    # Migrad flags
    "valid",
    "converged",
    "accurate_cov",
    "fmin_edm",

    # Migrad inputs
    "tolerance",
    "strategy",
    "ncall",

    # time
    "timestamp",
    "hash"
)

ARRAY_KEYS = (
    "fit_range",         # (n_pts,) timeslices fit
    "y_avg",             # (n_pts,) data central values
    "y_err",             # (n_pts,) data errors
    "cov",               # (n_pts, n_pts)
    "cov_inv",           # (n_pts, n_pts)
    "residuals",         # (n_pts,) y_data - y_model at central values
    "valid_res",         # (n_res,) bool — convergence per resample
    "t_dense",           # (n_dense,) fine grid for plotting
    "model_eval_res",    # (n_res, n_dense) fit evaluated per resample on fine grid
)

DICT_KEYS = (
    "p_cen",         # {name: float} — physical units
    "p_err",         # {name: float} — jackknife/bootstrap error
    "p_err_hesse",   # {name: float} — Hesse error (cross-check)
    "p_res",         # {name: (n_res,) array} — resample values, physical units
    "p_start",       # {name: float} — starting values (provenance)
    "p_limits",      # {name: (lo, hi)} — bounds (provenance)
)


def collect_quality(
    chi2:   float,
    n_par:  int,
    n_data: int,
) -> dict:

    n_dof = n_data - n_par

    return {
        "n_dof":    n_dof,
        "chi2_red": chi2 / n_dof,
        "pvalue":   float(chi2_dist.sf(chi2, n_dof)),
        "aic":      chi2 + 2 * n_par,
        "bic":      chi2 + n_par * np.log(n_data),
    }

def compute_residuals(y_data, y_model):
    """Return residuals and pulls (uncorrelated)."""
    residuals = y_data - y_model
    return residuals