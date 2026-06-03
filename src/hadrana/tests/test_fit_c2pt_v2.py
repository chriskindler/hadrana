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
    estimate_initial_timeslices,
    estimate_c2pt_minimum_timeslice,
    estimate_maximum_timeslice,
)

from hadrana.loader import load_c2pt_per_nsquare
from hadrana.runs.serialise import (
    Run,
    fit_results,
    make_fit_hash,
    make_c2pt_fit_id,
    save_fit
) 

def print_fit_summary(
    label: str,
    fit,
    path: Path = None,
    resample_method: str = "jackknife",
    accepted: bool | None = None,
    chi2dof_max: float = 2.0,
):
    c       = fit.central
    nres    = len(fit.resample) if fit.resample else 0
    nvalid  = sum(r.valid for r in fit.resample) if fit.resample else 0
    chi2dof = c.chi2 / c.ndof

    if accepted is None:
        accepted = c.judge_fit and chi2dof <= chi2dof_max

    if not c.judge_fit:
        status = "FAIL(nonfinite)"
    elif accepted is False:
        status = "REJECT(chi2)"
    elif c.edm > 1e-3 * c.tolerance:
        status = "ok*(edm)"
    else:
        status = "ok"

    # resample errors (matches fit_results convention)
    if fit.resample:
        factor = np.sqrt(nres - 1) if resample_method == "jackknife" else 1.0
        ddof   = 0 if resample_method == "jackknife" else 1
        params_err = {
            k: float(factor * np.std([r.params_est[k] for r in fit.resample], ddof=ddof))
            for k in c.params_est
        }
    else:
        params_err = c.params_err_hess

    bar = "-" * 64
    print(bar)
    print(f"  {label:<14s} [{status}]")
    print(f"    chi2/ndof = {chi2dof:7.3f}   (chi2 = {c.chi2:8.3f}, ndof = {c.ndof:3d})")
    print(f"    p-value   = {c.pvalue:7.3f}")
    print(f"    edm       = {c.edm:7.1e}   (tol = {c.tolerance:.1e})")
    print(f"    AIC/AICc  = {c.aic:7.2f} / {c.aicc:7.2f}")
    print(f"    resamples = {nvalid:4d}/{nres:<4d} valid")

    print(f"    {'param':>6s}  {'estimate':>14s}  {'jkn err':>12s}  {'hesse err':>12s}")
    for k, v in c.params_est.items():
        e_res   = params_err[k]
        e_hesse = c.params_err_hess[k]
        print(f"    {k:>6s}  {v:14.6e}  {e_res:12.3e}  {e_hesse:12.3e}   {gv.gvar(v, e_res)}")

    if path is not None:
        print(f"    saved: {path.name}")
    print(bar)

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

    # TODO: Find a way to make start params dict more flexible
    def exclude_keys(d, keys):
        return {x: d[x] for x in d if x not in keys}

    def accept_two_state(fit, chi2dof_max=2.0):
        c = fit.central
        return c.judge_fit and (c.chi2 / c.ndof) <= chi2dof_max

    ensemble        = "D251"
    observable      = "c2pt"
    resample_type   = "jackknife"
    bin_size        = 2
    momentum_shells = [0, 1, 2, 3, 4, 5, 6, 8]
    correlation_type = "correlated"
    snr_threshold   = 5

    """
        Construct dynamic array depending on lattice spacing of ensmeble data
        t_zero_max = t_phys // a 
    """

    c2pt_per_nsquare = load_c2pt_per_nsquare(ensemble, momentum_shells, bin_size)

    run = Run.initialise(
        base_path  = Path("/home/ck/phd/results"),
        ensemble   = ensemble,
        observable = observable,
        label      = "test-run23",
    )

    run_dir = run.create(exist_ok=True)
    hash_length = 6

    base_spec = {
        "ensemble":                 "D251",
        "resample_type":            "jackknife",
        "bin_size":                 bin_size,
        "correlation_type":         correlation_type,
        "execute_central_fit":      True,
        "starting_value_strategy": "data_driven",
        "tolerance":                1e-5,
        "strategy":                 2,
        "ncall":                    10000,
    }

    for nsquare in momentum_shells:
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

        initial_timeslices = estimate_initial_timeslices(
            a_fm   = 0.064,
            t_max  = t_max,
            t_phys = 0.45,
            t_start = 2,
            min_window = 4,
        )

        # ----------------------------------------------------------
        # PHASE 1: sweep t_zero, collect candidate t_min values
        # ----------------------------------------------------------
        candidate_tmins = []

        for t_zero in initial_timeslices:
            print(f"\n=== nsquare={nsquare}, t_zero={t_zero}, t_max={t_max} ===")

            params_start = estimate_c2pt_starting_values(
                c2pt_jkn = y,
                t_start  = t_zero,
                t_final  = t_max,
            )
            params_limit = {
                "A0":  (0, None),
                "E0":  (0, None),
                "A1":  (0, None),
                "dE1": (0, None),
            }

            fit_spec = {
                **base_spec,
                "correlation_type":     "correlated",
                "nsquare":              nsquare,
                "execute_resample_fit": True,
                "model_id":             "two-state-exp",
                "t_start":              t_zero,
                "t_final":              t_max,
                "bin_size":             bin_size,
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
                ncall                = fit_spec["ncall"],
            )

            fit_hash = make_fit_hash(fit_spec, hash_length)
            fit_id   = make_c2pt_fit_id(fit_spec, hash_length)
            result = fit_results(
                run_id          = run.run_id,
                fit_id          = fit_id,
                fit_hash        = fit_hash,
                fit_spec        = fit_spec,
                central         = two_state_fit.central,
                resample        = two_state_fit.resample,
                resample_method = fit_spec["resample_type"],
            )
            save_fit(result, run_dir, fit_id)
            acc = accept_two_state(two_state_fit)
            print_fit_summary("two-state", two_state_fit, accepted=acc)
            if not acc:
                print("two-state NOT accepted (chi2/ndof), skip")
                continue

            t_min = estimate_c2pt_minimum_timeslice(
                c2pt_jkn_err      = y_err,
                initial_timeslice = t_zero,
                maximum_timeslice = t_max,
                A1_est            = two_state_fit.central.params_est["A1"],
                E0_est            = two_state_fit.central.params_est["E0"],
                dE1_est           = two_state_fit.central.params_est["dE1"],
            )
            print(f"candidate t_min (t_zero={t_zero}): {t_min}")
            if t_min is not None:
                candidate_tmins.append(t_min)

        # ----------------------------------------------------------
        # COLLAPSE: largest t_min over all accepted two-state fits
        # ----------------------------------------------------------
        if not candidate_tmins:
            print(f"nsquare={nsquare}: no accepted two-state fit, skip one-state")
            continue

        t_min = max(candidate_tmins)
        print(f"nsquare={nsquare}: collapsed t_min = {t_min} "
            f"from candidates {sorted(candidate_tmins)}")

        # ----------------------------------------------------------
        # PHASE 2: single uncorrelated one-state production fit
        # ----------------------------------------------------------
        params_start = estimate_c2pt_starting_values(
            c2pt_jkn = y,
            t_start  = t_min,
            t_final  = t_max,
        )
        params_start = exclude_keys(params_start, {"A1", "dE1"})
        params_limit = {
            "A0": (0, None),
            "E0": (0, None),
        }

        fit_spec = {
            **base_spec,
            "correlation_type":     correlation_type,
            "nsquare":              nsquare,
            "execute_resample_fit": True,
            "model_id":             "one-state-exp",
            "t_start":              t_min,
            "t_final":              t_max,
            "bin_size":             bin_size,
            "params_start":         params_start,
            "params_limit":         params_limit,
        }

        one_state_fit = run_c2pt_fits(
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
            ncall                = fit_spec["ncall"],
        )

        fit_hash = make_fit_hash(fit_spec, hash_length)
        fit_id   = make_c2pt_fit_id(fit_spec, hash_length)
        result = fit_results(
            run_id          = run.run_id,
            fit_id          = fit_id,
            fit_hash        = fit_hash,
            fit_spec        = fit_spec,
            central         = one_state_fit.central,
            resample        = one_state_fit.resample,
            resample_method = fit_spec["resample_type"],
        )
        path_one = save_fit(result, run_dir, fit_id)
        print(f"binsize={bin_size}")
        print(f"FIT RANGE ONE-STATE = [{t_min}, {t_max}]")
        print_fit_summary("one-state", one_state_fit, path_one)