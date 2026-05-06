import logging
from pathlib import Path

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