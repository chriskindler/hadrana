# hadrana/fits/results.py
from __future__ import annotations
import numpy as np

from datetime import datetime, timezone
from typing import Literal, Optional
from scipy.stats import chi2 as chi2_dist

from chigrad.fit import FitResult

# The fit-result dict has three categories of values, each goes to a different
# place during serialization:
#   - SCALAR_KEYS  → manifest (parquet) and HDF5 attributes
#   - ARRAY_KEYS   → HDF5 datasets (too bulky for parquet)
#   - DICT_KEYS    → flattened into manifest columns; saved as nested groups in HDF5

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
    "params_res",         # {name: (n_res,) array} - resample values
    "params_cen",         # {name: float} - central values
    "params_err",         # {name: float} - jackknife/bootstrap error
    "params_err_hesse",   # {name: float} - Hesse error from migrad (cross-check with jackknife/bootstrap error)
    "params_start",       # {name: float} - starting values
    "params_limits",      # {name: (lo, hi)} - bounds
)

def fit_collect(
    t:               np.ndarray,
    y:               np.ndarray,
    central:         FitResult,
    resample:        list[FitResult] | None,
    resample_method: Literal["jackknife", "bootstrap"],
    fit_id:          str,
    fit_spec:        dict,
) -> dict:
    """Convert a (central, resamples) pair into the result-dict schema."""
    result = {
        "fit_id": fit_id,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "model_id": fit_spec["model_id"],
        "resample_method": resample_method,
        "correlation_type": fit_spec["correlation_type"],
        "bin_size": fit_spec["bin_size"],
        "t_min": fit_spec["t_min"],
        "t_max": fit_spec["t_max"],
        # counts
        "n_data": len(y[0]) if y.ndim == 2 else len(y),
        "n_par": central.npar,
        "n_dof": central.ndof,
        "n_res": len(resample) if resample else 0,
        "n_valid_res": sum(r.valid for r in resample) if resample else 0,
        # quality
        "chi2": central.chi2,
        "chi2_red": central.chi2_red,
        "pvalue": central.pvalue,
        "aic": central.aic,
        "bic": central.chi2 + central.npar * np.log(len(y[0]) if y.ndim == 2 else len(y)),
        # convergence
        "valid": central.valid,
        "converged": central.converged,
        "accurate_cov": central.accurate_covariance,
        "fmin_edm": central.fmin.edm,
        # optimizer settings (provenance)
        "tolerance": fit_spec.get("tolerance"),
        "strategy": fit_spec.get("strategy"),
        "ncall": fit_spec.get("ncall"),
        # arrays
        "fit_range": t,
        # dicts
        "p_cen": dict(central.values),
        "p_err_hesse": dict(central.errors),
        "p_limits": fit_spec.get("limits", {}),
    }

    if resample is not None:
        # per-parameter resample arrays
        param_names = list(central.values.keys())
        p_res = {name: np.array([r.values[name] for r in resample]) for name in param_names}
        result["p_res"] = p_res

        # jackknife/bootstrap errors
        if resample_method == "jackknife":
            n = len(resample)
            p_err = {name: float(np.sqrt(n - 1) * np.std(p_res[name], ddof=0)) for name in param_names}
        else: # bootstrap
            p_err = {name: float(np.std(p_res[name], ddof=1)) for name in param_names}
        result["p_err"] = p_err

        result["valid_res"] = np.array([r.valid for r in resample])

    return result