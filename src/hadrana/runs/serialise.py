# hadrana/runs/serialise.py
from __future__ import annotations
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal, Optional
import h5py
import hashlib
import json
import logging
import numpy as np
import pandas as pd

from chigrad.result import FitResult, FitRunResult

log = logging.getLogger(__name__)

Observable = Literal["c2pt", "c2pt_ratio", "ratio", "formfactor"]

class Run:
    def __init__(self, base_path: Path, ensemble: str, observable: Observable, run_id: str):
        self.base_path  = Path(base_path)
        self.ensemble   = ensemble
        self.observable = observable
        self.run_id     = run_id
        self.path       = self.base_path / observable / ensemble / run_id

    @classmethod
    def initialise(cls, base_path: Path, ensemble: str, observable: Observable, label: str) -> Run:
        run_id = f"{ensemble}-{observable}-fits-{label}"
        return cls(base_path, ensemble, observable, run_id)
    
    def create(self, exist_ok: bool = False) -> Path:
        self.path.mkdir(parents=True, exist_ok=exist_ok)
        return self.path

def fit_results(
    run_id:          str,
    fit_id:          str,
    fit_hash:        str,
    fit_spec:        dict,
    central:         FitResult,
    resample:        Optional[list[FitResult]],
    resample_method: Literal["jackknife", "bootstrap"],
) -> dict:

    y_fit_err_ext = None
    valid_res     = None
    params_err    = {}
    params_res    = {}

    if resample:
        nres   = len(resample)
        keys   = list(central.params_est.keys())
        factor = np.sqrt(nres - 1) if resample_method == "jackknife" else 1.0
        ddof   = 0 if resample_method == "jackknife" else 1

        y_fit_res_ext = np.array([r.y_fit_ext for r in resample])
        y_fit_err_ext = factor * np.std(y_fit_res_ext, axis=0, ddof=ddof)

        params_res = {k: np.array([r.params_est[k] for r in resample]) for k in keys}
        params_err = {k: float(factor * np.std(params_res[k], ddof=ddof)) for k in keys}

        valid_res  = np.array([r.valid for r in resample])

    result = {
        # identification
        "run_id":           run_id,
        "fit_id":           fit_id,
        "fit_hash":         fit_hash,
        "timestamp":        datetime.now(timezone.utc).strftime("%Y-%m-%d-UTC-T%H-%M-%S"),
        "ensemble":         fit_spec["ensemble"],
        "model_id":         fit_spec["model_id"],
        "resample_method":  resample_method,
        "correlation_type": fit_spec["correlation_type"],
        "nsquare":          fit_spec["nsquare"],
        "bin_size":         fit_spec["bin_size"],
        # counts
        "ndata":            central.ndata,
        "npar":             central.npar,
        "ndof":             central.ndof,
        "nres":             len(resample) if resample else 0,
        "nvalid":           sum(r.valid for r in resample) if resample else 0,
        # quality
        "chi2":             central.chi2,
        "chi2dof":          central.chi2 / central.ndof if central.ndof > 0 else np.nan,
        "pvalue":           central.pvalue,
        "aic":              central.aic,
        "aicc":             central.aicc,
        # diagnostic
        "valid":            central.valid,
        "converged":        central.converged,
        "accurate_cov":     central.accurate_cov,
        "fmin_edm":         central.edm,
        # meta
        "tolerance":        central.tolerance,
        "strategy":         central.strategy,
        "ncall":            central.ncall,
        # arrays
        "fit_range":        central.t,
        "fit_range_min":    central.t[0],
        "fit_range_max":    central.t[-1],
        "fit_range_ext":    central.t_ext,
        "residuals":        central.residuals,
        "y_fit":            central.y_fit,
        "y_fit_cen_ext":    central.y_fit_ext,
        "y_fit_err_ext":    y_fit_err_ext,
        "valid_res":        valid_res,
        # dicts
        "params_start":     fit_spec.get("params_start", {}) or {},
        "params_limit":     fit_spec.get("params_limit", {}) or {},
        "params_cen":       dict(central.params_est),
        "params_err":       params_err,
        "params_err_hesse": dict(central.params_err_hess),
        "params_res":       params_res,
    }

    return result
    
def make_fit_hash(fit_spec: dict, hash_length: int) -> str:
    """6-char hash of the full spec dict."""
    canonical = json.dumps(fit_spec, sort_keys=True, default=str)
    return hashlib.sha1(canonical.encode()).hexdigest()[:hash_length]

def make_c2pt_fit_id(fit_spec: dict, hash_length: int) -> str:
    """spec must contain: ensemble, t_zero, t_min, t_max, model_id,
    correlation_type, plus any optimizer settings.
    
    include hash in fit identifier as a safety net 
    """

    model_id = fit_spec["model_id"]
    elements = [
        fit_spec["ensemble"], "c2pt",
        f"nsquare{fit_spec["nsquare"]:02d}",
        f"tmin{fit_spec['t_start']:02d}",
        f"tmax{fit_spec['t_final']:02d}",
        model_id, fit_spec["correlation_type"]
    ]

    return "-".join(elements)

def make_c2pt_ratio_fit_id():
    pass

def make_ratio_fit_id():
    pass

def make_formfactor_fit_it():
    pass

# ---------------------------------------------------------------------------
# FIT I/O
# ---------------------------------------------------------------------------

SCALAR_KEYS = (
    "run_id", "fit_id", "fit_hash", "timestamp", "ensemble", "model_id",
    "resample_method", "correlation_type", "nsquare", "bin_size",
    "ndata", "npar", "ndof", "nres", "nvalid",
    "chi2", "chi2dof", "pvalue", "aic", "aicc",
    "valid", "converged", "accurate_cov", "fmin_edm",
    "tolerance", "strategy", "ncall",
    "fit_range_min", "fit_range_max",
)
ARRAY_KEYS = (
    "fit_range", "fit_range_ext", "residuals",
    "y_fit", "y_fit_cen_ext", "y_fit_err_ext", "valid_res",
)
DICT_KEYS = (
    "params_start", "params_limit",
    "params_cen", "params_err_hesse",
    "params_res", "params_err",
)

def save_fit(result: dict, run_dir: Path, fit_id: str) -> Path:
    """Write one fit's HDF5 artifact, sorted into /scalars, /arrays, /dicts."""
    path = run_dir / f"{fit_id}.h5"
    with h5py.File(path, "w") as f:

        # SCALARS
        scalars = f.create_group("scalars")
        for key in SCALAR_KEYS:
            if key not in result:
                continue
            v = result[key]
            if v is None:
                v = ""
            if isinstance(v, str):
                scalars.create_dataset(key, data=v, dtype=h5py.string_dtype())
            else:
                scalars.create_dataset(key, data=v)

        # ARRAYS
        arrays = f.create_group("arrays")
        for key in ARRAY_KEYS:
            if key not in result or result[key] is None:
                continue
            arrays.create_dataset(key, data=result[key], compression="gzip")

        # DICTS
        dicts = f.create_group("dicts")
        for key in DICT_KEYS:
            if key not in result:
                continue
            grp = dicts.create_group(key)
            for name, value in result[key].items():
                if isinstance(value, np.ndarray):
                    grp.create_dataset(name, data=value, compression="gzip")
                elif isinstance(value, tuple):  # bounds (lo, hi)
                    lo, hi = value
                    lo = -np.inf if lo is None else float(lo)
                    hi =  np.inf if hi is None else float(hi)
                    grp.attrs[name] = np.array([lo, hi], dtype=float)
                else:
                    grp.attrs[name] = value if value is not None else ""

    log.debug(f"saved fit: {path}")
    return path

# ---------------------------------------------------------------------------
# MANIFEST
# ---------------------------------------------------------------------------

def manifest_row(result: dict) -> dict:
    """Flatten a fit result dict into a flat row for the parquet manifest."""
    row = {k: result[k] for k in SCALAR_KEYS if k in result}
    # Flatten p_cen and p_err into per-parameter columns
    for name, value in result.get("p_cen", {}).items():
        row[f"{name}_cen"] = value
    for name, value in result.get("p_err", {}).items():
        row[f"{name}_err"] = value
    return row

def write_manifest(rows: list[dict], run_dir: Path) -> Path:
    out = run_dir / "fits.parquet"
    pd.DataFrame(rows).to_parquet(out, index=False)
    return out

if __name__ == "__main__":
    base_config = {
        "ensemble": "D251",
        "nsquare": 1,
        "bin_size": 1,
        "t_min": 13,
        "t_max": 44,
        "start_params": None,
        "model_id": "one-state-exp",
        "correlation_type": "uncorrelated",
        "tolerance": 1e-4,
        "strategy": 1,
        "ncall": 10000,
    }

    run_id = "run03"

    ensemble   = "D251"
    observable = "c2pt" 
    base_path  = Path("/home/ck/phd/results")

    hash_length: int = 12

    # run     = Run.initialise(base_path, ensemble, observable, run_id)
    # run_dir = run.create(exist_ok=False)

    # print(run_dir)

    model_id = "one-state-exp"
    initial_timeslices = [3, 4, 5, 6, 7, 8]

    """
    def make_c2pt_fit_id(
    ensemble:    str,
    fit_config:  dict,
    hash_length: int,
    t_zero:      Optional[int] = None,
)
    """

    for t_zero in initial_timeslices:
        fit_config = {**base_config, "t_zero": t_zero}
        fit_id = make_c2pt_fit_id(fit_config)
        print(fit_id)