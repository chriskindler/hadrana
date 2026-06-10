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

from hadrana.loader import load_c2pt_per_nsquare, load_rwfs
from hadrana.runs.serialise import (
    Run,
    fit_results,
    make_fit_hash,
    make_c2pt_fit_id,
    save_fit
) 

from hadrana.statistics import maximum_binsize

def print_fit_summary(
    label: str,
    fit,
    path: Path,
    resample_method: str = "jackknife",
    accepted: bool | None = None,
    chi2dof_max: float = 2.0,
):
    c       = fit.central
    nres    = len(fit.resample) if fit.resample else 0
    nvalid  = sum(r.valid for r in fit.resample) if fit.resample else 0

    # guard against degenerate fit windows (ndof <= 0)
    ndof_ok = c.ndof > 0
    chi2dof = (c.chi2 / c.ndof) if ndof_ok else float("inf")

    if accepted is None:
        accepted = ndof_ok and c.judge_fit and chi2dof <= chi2dof_max

    if not ndof_ok:
        status = "FAIL(ndof<=0)"
    elif not c.judge_fit:
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
        if c.ndof <= 0:
            return False
        return c.judge_fit and (c.chi2 / c.ndof) <= chi2dof_max

    def jackknife_err(fit, key, resample_method="jackknife"):
        """jackknife sigma of a single parameter's resample point estimates."""
        if not fit.resample:
            return float(fit.central.params_err_hess[key])
        nres   = len(fit.resample)
        factor = np.sqrt(nres - 1) if resample_method == "jackknife" else 1.0
        ddof   = 0 if resample_method == "jackknife" else 1
        return float(factor * np.std([r.params_est[key] for r in fit.resample], ddof=ddof))

    ensemble        = "A654"
    observable      = "c2pt"
    resample_type   = "jackknife"
    momentum_shells = [2, 3, 4, 5, 6, 8]
    snr_threshold   = 10
    a_fm            = 0.098
    t_phys          = 0.3
    t_start         = 2

    # binsize used to SELECT and LOCK the fit range (small => many bins =>
    # well-conditioned correlated covariance => chi2/ndof <= 2.0 is meaningful)
    selection_binsize = 2

    rwfs, rwfs_path = load_rwfs(ensemble)
    ncfg = rwfs.shape[0]

    binsizes = [1, 2, 4, 8, 16, 32, 64, maximum_binsize(ncfg)]
    # de-duplicate (maximum_binsize may collide with a literal) while preserving order
    binsizes = sorted(set(binsizes))

    print(f"BINSIZES = {binsizes}")
    print(f"maximum binsize = {np.max(binsizes)}")
    print(f"selection binsize = {selection_binsize}")

    # one run dir for the whole study
    run = Run.initialise(
        base_path  = Path("/home/ck/phd/results"),
        ensemble   = ensemble,
        observable = observable,
        label      = "binsize-run01",
    )
    run_dir     = run.create(exist_ok=True)
    hash_length = 6

    base_spec = {
        "ensemble":                ensemble,
        "resample_type":           "jackknife",
        "execute_central_fit":     True,
        "starting_value_strategy": "data_driven",
        "tolerance":               1e-5,
        "strategy":                2,
        "ncall":                   10000,
    }

    # ==============================================================
    # SELECTION PHASE (run ONCE at selection_binsize)
    #   correlated two-state sweep over t_zero, chi2/ndof <= 2.0 guard,
    #   collapse to largest accepted t_min; lock (t_min, t_max) per shell.
    # ==============================================================
    print("\n" + "=" * 64)
    print(f"SELECTION PHASE @ binsize={selection_binsize}")
    print("=" * 64)

    c2pt_sel = load_c2pt_per_nsquare(ensemble, momentum_shells, selection_binsize)

    # nsquare -> (t_min, t_max); shells with no accepted two-state are absent
    locked_ranges: dict[int, tuple[int, int]] = {}

    for nsquare in momentum_shells:
        y     = c2pt_sel[nsquare]
        y_cen = np.mean(y, axis=0)
        y_err = compute_sdev(y, resample_type)

        # MAXIMUM TIMESLICE (locked here at the selection binsize)
        t_max = estimate_maximum_timeslice(
            c2pt_jkn_avg              = y_cen,
            c2pt_jkn_err              = y_err,
            signal_to_noise_threshold = snr_threshold,
        )

        print(f"tmax = {t_max}")

        initial_timeslices = estimate_initial_timeslices(
            a_fm    = a_fm,
            t_phys = t_phys,
            t_max  = t_max,
        )

        candidate_tmins = []

        for t_zero in initial_timeslices:
            print(f"\n=== [select] nsquare={nsquare}, t_zero={t_zero}, t_max={t_max} ===")

            # need at least (nparams + 1) points for a meaningful two-state ndof
            if (t_max - t_zero + 1) <= 4:
                print(f"window [{t_zero},{t_max}] too short for two-state (<=4 pts), skip")
                continue

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
                "binsize":              selection_binsize,
                "correlation_type":     "correlated",
                "nsquare":              nsquare,
                "execute_resample_fit": True,
                "model_id":             "two-state-exp",
                "t_start":              t_zero,
                "t_final":              t_max,
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
            path_two = save_fit(result, run_dir, fit_id)
            acc = accept_two_state(two_state_fit)
            print_fit_summary("two-state", two_state_fit, path=path_two, accepted=acc)
            if not acc:
                print("two-state NOT accepted (chi2/ndof or ndof<=0), skip")
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

        if not candidate_tmins:
            print(f"nsquare={nsquare}: no accepted two-state fit, NOT locking range")
            continue

        t_min = max(candidate_tmins)
        locked_ranges[nsquare] = (int(t_min), int(t_max))
        print(f"nsquare={nsquare}: LOCKED range [{t_min}, {t_max}] "
              f"from candidates {sorted(candidate_tmins)}")

    print("\n" + "=" * 64)
    print("LOCKED RANGES")
    for ns in momentum_shells:
        if ns in locked_ranges:
            print(f"  nsquare={ns}: {locked_ranges[ns]}")
        else:
            print(f"  nsquare={ns}: (none - excluded from sweep)")
    print("=" * 64)

    # ==============================================================
    # BINSIZE SWEEP (uncorrelated one-state over the LOCKED ranges)
    #   jackknife sigma(E0) is invariant to correlated/uncorrelated;
    #   uncorrelated avoids the singular-covariance blowup at large binsize.
    # ==============================================================
    # nsquare -> list of (binsize, E0, sigmaE0) for the saturation curve
    saturation: dict[int, list[tuple[int, float, float]]] = {
        ns: [] for ns in locked_ranges
    }

    for s in binsizes:
        print("\n" + "#" * 64)
        print(f"BINSIZE SWEEP @ binsize={s}")
        print("#" * 64)

        c2pt_per_nsquare = load_c2pt_per_nsquare(ensemble, momentum_shells, s)

        for nsquare in momentum_shells:
            if nsquare not in locked_ranges:
                continue

            t_min, t_max = locked_ranges[nsquare]
            y = c2pt_per_nsquare[nsquare]

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
                "binsize":              s,
                "correlation_type":     "uncorrelated",
                "nsquare":              nsquare,
                "execute_resample_fit": True,
                "model_id":             "one-state-exp",
                "t_start":              t_min,
                "t_final":              t_max,
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

            E0      = float(one_state_fit.central.params_est["E0"])
            sigmaE0 = jackknife_err(one_state_fit, "E0", resample_type)
            saturation[nsquare].append((s, E0, sigmaE0))

            print(f"\nbinsize={s}  nsquare={nsquare}")
            print(f"FIT RANGE ONE-STATE = [{t_min}, {t_max}]")
            print_fit_summary("one-state", one_state_fit, path_one)

    # ==============================================================
    # SATURATION SUMMARY: sigma(E0) vs binsize, per momentum shell
    # ==============================================================
    print("\n" + "=" * 64)
    print("SATURATION CURVE: sigma(E0) vs binsize")
    print("=" * 64)
    for nsquare in momentum_shells:
        if nsquare not in saturation:
            continue
        t_min, t_max = locked_ranges[nsquare]
        print(f"\nnsquare={nsquare}  range=[{t_min},{t_max}]")
        print(f"  {'binsize':>8s}  {'E0':>14s}  {'sigma(E0)':>14s}")
        for s, E0, sE0 in saturation[nsquare]:
            print(f"  {s:8d}  {E0:14.6e}  {sE0:14.6e}")
