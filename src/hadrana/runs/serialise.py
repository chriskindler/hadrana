# hadrana/runs/serialise.py
from __future__ import annotations
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal
import h5py
import hashlib
import json
import logging
import numpy as np
import pandas as pd

from hadrana.fits.results import SCALAR_KEYS, ARRAY_KEYS, DICT_KEYS

log = logging.getLogger(__name__)

"""
base path: /home/ck/phd/code/nucleon_analysis/results/c2pt_fits/nmax0{nmax}
directory name = run_dir
{timestamp}-{ensemble}-{observable}-nmax0{nmax}-{label}

/home/ck/phd/code/nucleon_analysis/results/c2pt_fits/mom01
    2026-05-04-T14-20-D251-c2pt-mom01-bin01-snr05-one-state/
        D251-c2pt-mom01-tzero04-tmin16-tmax44-one-state-correlated.h5
        D251-c2pt-mom01-tzero03-tmin18-tmax44-one-state-uncorrelated.h5

    2026-05-04-T14-20-D251-c2pt-mom01-bin01-snr05-two-state/mom01
        2026-05-04-T14-20-D251-c2pt-mom01-tzero04-tmin16-tmax44-two-state-correlated.h5
        2026-05-04-T14-20-D251-c2pt-mom01-tzero03-tmin18-tmax44-two-state-uncorrelated.h5

/home/ck/phd/code/nucleon_analysis/results/c2pt_ratio_fits/mom01
    2026-05-04-T14-20-D251-c2pt-ratio-nmax01-constant/
        D251-c2pt-ratio-mom01-tzero04-tmin16-tmax44-constant-correlated-a3f2b1.h5
        D251-c2pt-ratio-mom01-tzero04-tmin16-tmax44-constant-uncorrelated-a3f2b1.h5

    2026-05-04-T14-20-D251-c2pt-ratio-nmax01-excited-overlap-included/
        D251-c2pt-ratio-mom01-tzero04-tmin16-tmax44-excited-overlap-included-correlated-a3f2b1.h5
        D251-c2pt-ratio-mom01-tzero04-tmin16-tmax44-excited-overlap-included-uncorrelated-a3f2b1.h5

    2026-05-04-T14-20-D251-c2pt-ratio-nmax01-excited-overlap-excluded/
        D251-c2pt-ratio-mom01-tzero04-tmin16-tmax44-excited-overlap-excluded-correlated-a3f2b1.h5
        D251-c2pt-ratio-mom01-tzero04-tmin16-tmax44-excited-overlap-excluded-uncorrelated-a3f2b1.h5
"""


# ---------------------------------------------------------------------------
# IDS
# ---------------------------------------------------------------------------

def make_fit_hash(spec: dict, length: int = 6) -> str:
    """6-char hash of the full spec dict."""
    canonical = json.dumps(spec, sort_keys=True, default=str)
    return hashlib.sha1(canonical.encode()).hexdigest()[:length]

def make_fit_id(spec: dict) -> tuple[str, str]:
    """Composite human-readable ID + 6-char hash suffix."""
    canonical = json.dumps(spec, sort_keys=True, default=str)
    fit_hash = hashlib.sha1(canonical.encode()).hexdigest()[:6]
    parts = [
        spec["ensemble"],
        f"nmax{spec['nmax']:02d}",
        f"tmin{spec['t_min']:02d}",
        f"tmax{spec['t_max']:02d}",
        spec["model_id"],
        fit_hash,
    ]
    return "-".join(parts), fit_hash

def make_run_id(ensemble: str, observable: str, label: str) -> str:
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%MZ")
    return f"{timestamp}_{ensemble}_{observable}_{label}"

def make_run_directory(
    analysis_root_path: Path,
    ensemble:           str,
    observable:         Literal["c2pt", "c2pt_ratio", "ratio", "formfactor"],
    label:              str,
) -> Path:
    """Create a fresh run directory for one campaign."""
    run_id = make_run_id(ensemble, observable, label)
    run_dir = analysis_root_path / observable / "runs" / run_id
    run_dir.mkdir(parents=True, exist_ok=False)
    (run_dir / "fits").mkdir()
    (run_dir / "inputs").mkdir()
    (run_dir / "plots").mkdir()
    return run_dir

# ---------------------------------------------------------------------------
# FIT SAVE / LOAD
# ---------------------------------------------------------------------------

def save_fit(result: dict, run_dir: Path, fit_id: str) -> Path:
    """Write one fit's HDF5 artifact."""
    path = run_dir / "fits" / f"{fit_id}.h5"
    with h5py.File(path, "w") as f:
        # Scalars as root attrs
        for key in SCALAR_KEYS:
            if key in result:
                v = result[key]
                f.attrs[key] = "" if v is None else v

        # Arrays as datasets
        for key in ARRAY_KEYS:
            if key in result:
                f.create_dataset(key, data=result[key], compression="gzip")

        # Dicts as groups with attrs (or datasets for resample arrays)
        for key in DICT_KEYS:
            if key not in result:
                continue
            grp = f.create_group(key)
            for name, value in result[key].items():
                if isinstance(value, np.ndarray):
                    grp.create_dataset(name, data=value, compression="gzip")
                elif isinstance(value, tuple):  # e.g., limits (lo, hi)
                    grp.attrs[name] = list(value)
                else:
                    grp.attrs[name] = value if value is not None else ""

    log.debug(f"saved fit: {path}")
    return path


def load_fit(fit_id: str, run_dir: Path) -> dict:
    """Load a fit's HDF5 artifact back to a dict."""
    path = run_dir / "fits" / f"{fit_id}.h5"
    with h5py.File(path, "r") as f:
        result = {key: f.attrs[key] for key in f.attrs}
        for key in ARRAY_KEYS:
            if key in f:
                result[key] = f[key][...]
        for key in DICT_KEYS:
            if key in f:
                grp = f[key]
                d = dict(grp.attrs)
                for name in grp:
                    d[name] = grp[name][...]
                result[key] = d
    return result


# ---------------------------------------------------------------------------
# MANIFEST
# ---------------------------------------------------------------------------

# def manifest_row(result: dict) -> dict:
#     """Flatten a fit result dict into a flat row for the parquet manifest."""
#     row = {k: result[k] for k in SCALAR_KEYS if k in result}
#     # Flatten p_cen and p_err into per-parameter columns
#     for name, value in result.get("p_cen", {}).items():
#         row[f"{name}_cen"] = value
#     for name, value in result.get("p_err", {}).items():
#         row[f"{name}_err"] = value
#     return row


# def write_manifest(rows: list[dict], run_dir: Path) -> Path:
#     out = run_dir / "fits.parquet"
#     pd.DataFrame(rows).to_parquet(out, index=False)
#     return out