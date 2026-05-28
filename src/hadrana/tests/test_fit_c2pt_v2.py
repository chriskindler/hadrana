# hadrana/tests/test_fit_c2pt_v2.py
import argparse
import gvar as gv
import h5py
import json
import pickle
import numpy as np

from iminuit.util import describe
from pathlib import Path
from scipy.stats import chi2 as chi2_dist
from typing import Callable, Literal, Optional

from chigrad.statistics import compute_cov, compute_cov_inv, compute_sdev
from chigrad.fit import execute_fits

from hadrana.fits.fit_c2pt.models_c2pt import (
    get_model,
    estimate_c2pt_starting_values,
)
from hadrana.fits.fit_c2pt.timeslice_criteria import (
    estimate_c2pt_minimum_timeslice,
    estimate_maximum_timeslice,
)

from hadrana.loader import load_c2pt_per_nsquare

from hadrana.runs.serialise import SCALAR_KEYS, save_fit

from hadrana.runs.serialise import (
    Run,
    fit_results,
    make_fit_hash,
    make_c2pt_fit_id,
    save_fit
) 

"""
1. Write c2pt fit wrapper
2. Write full perform_c2pt_fit main function
    PHASE 0:
        Initialise run
        Read c2pt base config via argparse
        Construct fit config dictionary

    PHASE 1:
    Load source-momentum averaged forward-limit c2pt jackknife data from disk
    Determine t_max with signal-to-noise threshold from c2pt jackknife data
        for t_zero in initial_timeslices:
            Estimate starting values from c2pt jackknife data
            Perform central two-state fit per nsquare in the window [t_zero, t_max]
                return A1, E0 and dE1
            Determine t_min from excited-state criterion 

    PHASE 2:
        for t_zero in initial_timeslices:
            Perform central and resample one-state fits per nsquare in the window [t_min, t_max]
                return {A0_res, E0_res}
            Compute central values and errors from {A0_res, E0_res}
                
    PHASE 3:
        Serialise fit results

"""

def run_c2pt_fits(
    t_start:              int,
    t_final:              int,
    y:                    np.ndarray, # (nres, ndata)

    correlation_type:     Literal["correlated", "uncorrelated"],
    resample_type:        Literal["bootstrap", "jackknife"],

    execute_central_fit:  bool,
    execute_resample_fit: bool,

    model_id:             str,
    params_start:         dict[str, float],
    params_limit:         Optional[dict[str, tuple[float | None, float | None]]],

    tolerance:            float,
    strategy:             int,
    ncall:                int,
):

    t     = np.arange(t_start, t_final + 1)
    t_ext = np.linspace(t_start, t_final, 500)

    y = y[:, t]

    cov_inv  = None
    sdev_inv = None

    if correlation_type == "correlated":
        cov     = compute_cov(y, resample_type)
        cov_inv = compute_cov_inv(cov)
    
    else:
        sdev     = compute_sdev(y, resample_type)
        sdev_inv = 1.0 / sdev

    model, _ = get_model(model_id)

    return execute_fits(
        t,
        t_ext,
        y,
        correlation_type,
        resample_type,
        cov_inv,
        sdev_inv,
        execute_central_fit,
        execute_resample_fit,
        model,
        params_start,
        params_limit,
        tolerance,
        strategy,
        ncall
    )

if __name__ == "__main__":
    # parser = argparse.ArgumentParser()
    # parser.add_argument("-f", "--filepath", help="json file path")
    # args = parser.parse_args()

    # with open(args.filepath, "r") as f:
    #     cfg: dict = json.load(f)

    print("SCALAR_KEYS:", SCALAR_KEYS)

    ensemble    = "D251"
    observable  = "c2pt"
    resample_type = "jackknife"
    nsquare_values = [0]
    initial_timeslices = [3, 4, 5, 6, 7, 8]
    snr_threshold  = 5
    t_zero = 5

    c2pt_per_nsquare = load_c2pt_per_nsquare(ensemble, nsquare_values)

    run = Run.initialise(
        base_path  = Path("/home/ck/phd/results"),
        ensemble   = ensemble,
        observable = observable,
        label      = "test-run10",
    )

    run_dir = run.create(exist_ok=True)
    hash_length = 6

    for nsquare in nsquare_values:
        
        #nmax         = 0
        y     = c2pt_per_nsquare[nsquare]
        nres  = y.shape[0]

        y_cen = np.mean(y, axis=0)
        y_err = compute_sdev(y, resample_type) 

        # MAXIMUM TIMESLICE
        t_max = estimate_maximum_timeslice(
            c2pt_jkn_avg              = y_cen,
            c2pt_jkn_err              = y_err,
            signal_to_noise_threshold = snr_threshold,
        )

        # FIT PARAMETERS
        params_start = estimate_c2pt_starting_values(
            c2pt_jkn = y,
            t_start  = t_zero,
            t_final  = t_max
        )

        params_limit = {
            "A0":  (0, None),
            "E0":  (0, None),
            "A1":  (0, None),
            "dE1": (0, None),
        }

        base_spec = {
            # general
            "ensemble":                 "D251",
            "resample_type":            "jackknife",
            "bin_size":                 1,
            "correlation_type":         "uncorrelated",
            "execute_central_fit":      True,
            "starting_value_strategy": "data_driven",
            "tolerance":                1e-4,
            "strategy":                 1,
            "ncall":                    10000,
        }

        fit_spec = {
            **base_spec,
            "nsquare":              nsquare,
            "execute_resample_fit": True,
            "model_id":             "two-state-exp",
            "t_start":              t_zero,
            "t_final":              t_max,
            "bin_size":             1,
            "params_start":         params_start,
            "params_limit":         params_limit,
        }

        two_state_fit = run_c2pt_fits(
            t_start              = fit_spec["t_start"],
            t_final              = fit_spec["t_final"],
            y                    = y,
            correlation_type     = fit_spec["correlation_type"],
            resample_type        = fit_spec["resample_type"],
            execute_central_fit  = fit_spec["execute_central_fit"],
            execute_resample_fit = fit_spec["execute_resample_fit"],
            model_id             = fit_spec["model_id"],
            params_start         = fit_spec["params_start"],
            params_limit         = fit_spec["params_limit"],
            tolerance            = fit_spec["tolerance"], 
            strategy             = fit_spec["strategy"],
            ncall                = fit_spec["ncall"]
        )

        print(two_state_fit.central.params_est)

        fit_hash = make_fit_hash(fit_spec, hash_length)
        fit_id   = make_c2pt_fit_id(fit_spec, hash_length)

        result = fit_results(
            run_id          = run.run_id,
            fit_id          = fit_id,
            fit_hash        = fit_hash,
            fit_spec        = fit_spec,
            central         = two_state_fit.central,
            resample        = two_state_fit.resample, # None — execute_resample_fit=False
            resample_method = fit_spec["resample_type"],
        )

        save_fit(result, run_dir, fit_id)

        A1_est  = two_state_fit.central.params_est["A1"]
        E0_est  = two_state_fit.central.params_est["E0"]
        dE1_est = two_state_fit.central.params_est["dE1"]
        

        t_min = estimate_c2pt_minimum_timeslice(
            c2pt_jkn_err      = y_err,
            initial_timeslice = t_zero,
            maximum_timeslice = t_max,
            A1_est            = A1_est,
            E0_est            = E0_est,
            dE1_est           = dE1_est,
        )


