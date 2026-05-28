# hadrana/runs/log.py
# DEFAULT
import logging
import sys
import socket
import getpass
import platform
import importlib.metadata
from pathlib import Path
from datetime import datetime, timezone

# CUSTOM
from hadrana.ensembles.ensemble_specs import ENSEMBLES

LOG_FORMAT = "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
LOG_DATEFMT = "%Y-%m-%dT%H:%M:%SZ"

def create_run_directory(ensemble: str, observable: str, run_id: str):
    """
    /home/ck/phd/results/c2pt_fits/D251
        2026-04-30T14-22-D251-c2pt-fits-run-id/
            config.json
            run_id.log
            fit_results/
                nmax00_tzero04_tmax44.h5
                    ├── /attrs                    # all manifest scalars duplicated as attrs
                    ├── /data
                    │   ├── c2pt                  (n_jkn, n_t)        # the input
                    │   ├── timeslices            (n_fit,)
                    │   ├── timeslices_extended   (n_density,)
                    ├── /fit
                    │   ├── A0_jkn                (n_jkn,)
                    │   ├── E0_jkn                (n_jkn,)
                    │   ├── valid_jkn             (n_jkn,) bool
                    │   ├── model_eval_jkn_extended  (n_jkn, n_density)
                    │   ├── cov                   (n_fit, n_fit)
                    │   └── cov_inv               (n_fit, n_fit)
                    └── /two_state_fit/...        # nested group with the seed fit
                nmax00_tzero05_tmax44.h5
                nmax00_tzero06_tmax44.h5
                nmax00_tzero07_tmax44.h5
        index.parquet
    """
    pass

def log_setup(observable: str, run_dir: Path):
    pass

def log_run_start():
    """ string for log file start"""
    pass

def log_ensemble_specs(ensemble: str):
    pass

def log_metadata(ensemble: str, observable: str):
    pass

def log_fit_diagnostic(ensemble: str, observable: str):
    pass

def log_summary(ensemble: str, observable: str):
    """log most important information and"""
    pass

def log_run_end(ensemble: str, observable: str):
    """ string for log file end, file export"""
    pass


"""
    // PHASE MARKERS

    PHASE A: RUN START

    PHASE 1: ENSEMBLE SPECS
        - id
        - trajectory
        - beta
        - lattice spacing
        - lattice dimension
        - boundary conditions in time
        - pion mass, kaon mass, LMpi, 
        - original config number

    PHASE 2: METADATA
        - run-id
        - timestamp
        - hostname
        - username
        - code version (cora, hadrana) (git hash)
        - import lib versions
        - path to fit config file
        - output fit config file for observable
        
    PHASE 3: DATA
        - input data filepath (+ checksums)
        - exceptional removal + indices

    PHASE 4: FIT DIAGNOSTIC
        - central value fit: summary
        - resample fit: only aggregate
        - table of results per nmax, t0 for c2pt fits 
        - fit summary

    PHASE 5: RUN SUMMARY
        - Attempted / Succeeded / Skipped / Failed
        - Total duration 
        - Path to file output

    PHASE Z: RUN END

"""

"""
DEBUG — things you'd only want when something went wrong:

Per-resample fit parameters
Intermediate matrix conditions, eigenvalue spectra
Detailed iteration counts inside Migrad

INFO — the narrative of the run; readable as a story:

Run start/end with config summary
Each fit being attempted with its specification
Each fit completing with key results
Phase transitions ("loading data", "starting fit campaign", "writing manifest")

WARNING — something unexpected but recoverable:

A fit didn't converge but campaign continues
Excited-state criterion never met (your existing print statement)
A resample fraction below threshold

ERROR — a fit/operation failed but the run continues:

Exception caught around a single fit
File missing for one nmax value

CRITICAL — the run cannot continue:

Config invalid, can't load raw data, etc.
"""

"""Logging configuration for hadrana analysis campaigns."""

import logging
from pathlib import Path
from typing import Optional

_FILE_FORMAT = (
    "%(asctime)s.%(msecs)03d | %(levelname)-7s | "
    "%(name)-24s | %(message)s"
)
_FILE_DATEFMT = "%Y-%m-%d %H:%M:%S"
_CONSOLE_FORMAT = "%(levelname)-7s | %(name)-24s | %(message)s"


def setup_logging(
    run_dir: Path,
    console_level: int = logging.INFO,
    file_level: int = logging.DEBUG,
    log_filename: str = "run.log",
) -> logging.Logger:
    """Configure the 'hadrana' logger for a campaign run.

    Writes everything at file_level+ to run_dir/run.log, and shows
    console_level+ on stderr. Idempotent — safe to call twice.
    """
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    file_handler = logging.FileHandler(run_dir / log_filename, mode="w")
    file_handler.setFormatter(logging.Formatter(_FILE_FORMAT, _FILE_DATEFMT))
    file_handler.setLevel(file_level)

    console_handler = logging.StreamHandler()
    console_handler.setFormatter(logging.Formatter(_CONSOLE_FORMAT))
    console_handler.setLevel(console_level)

    root = logging.getLogger("hadrana")
    root.setLevel(min(console_level, file_level))
    root.handlers.clear()
    root.addHandler(file_handler)
    root.addHandler(console_handler)
    root.propagate = False

    # also capture cora's logs into the same file
    cora_logger = logging.getLogger("cora")
    cora_logger.setLevel(file_level)
    cora_logger.handlers.clear()
    cora_logger.addHandler(file_handler)
    cora_logger.addHandler(console_handler)
    cora_logger.propagate = False

    return root


def get_logger(name: str) -> logging.Logger:
    """Convenience wrapper. Equivalent to logging.getLogger(name)."""
    return logging.getLogger(name)

if __name__ == "__main__":
    print("TEST LOG FILE GENERATION")
