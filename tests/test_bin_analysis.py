# hadrana/tests/test_bin_analysis.py
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

"""
    STRATEGY: D251
    binsizes = [1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20] 
    maxbin   = 20

    bin-run01
    momentum_shells = [0] later: [0, 1, 2, 3, 4, 5, 6, 8] 

    for n in momentum_shells:
        for s in binsizes:
            jkn = compute_jkn(n, s)
            avg = np.mean(jkn)
            err = np.std(jkn) * np.sqrt(ncfg - 1)

            tmax = determine_tmax(avg, err, snr_threshold)

"""

def print_fit_summary(
    label: str,
    fit,
    path: Path,
    resample_method: str = "jackknife",
    accepted: bool | None = None,
    chi2dof_max: float = 2.5,
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

    ensemble      = "A650"
    observable    = "c2pt"
    resample_type = "jackknife"

    momentum_shells = [0]
    t_phys          = 0.3
    t_start         = 2

    snr_threshold = 0
    a_fm = 0.0
    binsizes = []
    if ensemble == "D251":
        binsizes      = [1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
        snr_threshold = 10
        a_fm          = 0.064

    elif ensemble == "A650":
        binsizes      = [1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 30, 40, 50]
        snr_threshold = 5
        a_fm          = 0.098

    run = Run.initialise(
        base_path  = Path("/home/ck/phd/results"),
        ensemble   = ensemble,
        observable = observable,
        label      = "nucleon-mass-run01",
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

    # nsquare -> (t_min, t_max); shells with no accepted two-state are absent

    selection: dict[tuple[int, int], tuple[int, int]] = {}   # (nsquare, binsize) -> (t_min, t_max)
    records:   list[dict] = []                                   # one row per (n, s, t_zero)

    for n in momentum_shells:
        for s in binsizes:
            c2pt_per_nsquare = load_c2pt_per_nsquare(ensemble, momentum_shells, s)

            jkn = c2pt_per_nsquare[n]
            avg = np.mean(jkn, axis=0)
            err = compute_sdev(jkn, resample_type)

            # MAXIMUM TIMESLICE (locked here at the selection binsize)
            t_max = estimate_maximum_timeslice(
                c2pt_jkn_avg              = avg,
                c2pt_jkn_err              = err,
                signal_to_noise_threshold = snr_threshold,
            )

            print(f"tmax = {t_max}")

            initial_timeslices = estimate_initial_timeslices(
                a_fm   = a_fm,
                t_phys = t_phys,
                t_max  = t_max,
            )

            candidates = []

            for t_zero in initial_timeslices:
                print(f"\n [select] nsquare={n}, t_zero={t_zero}, t_max={t_max}")

                params_start = estimate_c2pt_starting_values(
                    c2pt_jkn = jkn,
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
                    "binsize":              s,
                    "correlation_type":     "correlated",
                    "nsquare":              n,
                    "execute_resample_fit": False,
                    "model_id":             "two-state-exp",
                    "t_start":              t_zero,
                    "t_final":              t_max,
                    "params_start":         params_start,
                    "params_limit":         params_limit,
                }

                two_state_fit = run_c2pt_fits(
                    t_start              = fit_spec["t_start"],
                    t_final              = fit_spec["t_final"],
                    y                    = jkn,
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

                c       = two_state_fit.central
                chi2dof = (c.chi2 / c.ndof) if c.ndof > 0 else float("inf")
                rec = {
                    "nsquare": n, "binsize": s, "t_zero": int(t_zero), "t_max": int(t_max),
                    "accepted": bool(acc), "chi2": float(c.chi2), "ndof": int(c.ndof),
                    "chi2dof": float(chi2dof), "pvalue": float(c.pvalue),
                    "t_min_candidate": None, "fit_id": fit_id,
                }

                if not acc:
                    print("two-state NOT accepted (chi2/ndof or ndof<=0), skip")
                    records.append(rec)
                    continue

                t_min_cand = estimate_c2pt_minimum_timeslice(
                    c2pt_jkn_err      = err,
                    initial_timeslice = t_zero,
                    maximum_timeslice = t_max,
                    A1_est            = two_state_fit.central.params_est["A1"],
                    E0_est            = two_state_fit.central.params_est["E0"],
                    dE1_est           = two_state_fit.central.params_est["dE1"],
                )
                print(f"candidate t_min (t_zero={t_zero}): {t_min_cand}")
                rec["t_min_candidate"] = None if t_min_cand is None else int(t_min_cand)
                records.append(rec)
                if t_min_cand is not None:
                    candidates.append(t_min_cand)

            if not candidates:
                print(f"nsquare={n}, binsize={s}: no accepted two-state fit, NOT locking range")
                continue

            t_min = max(candidates)
            selection[(n, s)] = (int(t_min), int(t_max))
            print(f"nsquare={n} binsize={s}: LOCKED [{t_min}, {t_max}] from {sorted(candidates)}")

    def summarise_selection(records, selection):
        by_bin: dict[tuple[int, int], list[dict]] = {}
        for r in records:
            by_bin.setdefault((r["nsquare"], r["binsize"]), []).append(r)

        print("\n" + "=" * 78)
        print("SELECTION SUMMARY")
        print("=" * 78)
        for (n, s) in sorted(by_bin):
            recs  = sorted(by_bin[(n, s)], key=lambda r: r["t_zero"])
            t_max = recs[0]["t_max"]
            lock  = selection.get((n, s))
            lock_s = f"[{lock[0]}, {lock[1]}]" if lock else "NOT LOCKED"
            print(f"\n nsquare={n}  binsize={s:>2}  t_max={t_max}  ->  locked {lock_s}")
            print(f"   {'t_zero':>6} {'acc':>4} {'chi2/ndof':>10} {'chi2':>9} {'ndof':>5} {'p':>6} {'t_min*':>7}")
            for r in recs:
                flag  = "Y" if r["accepted"] else "n"
                tmc   = r["t_min_candidate"]
                tmc_s = "-" if tmc is None else str(tmc)
                print(f"   {r['t_zero']:>6} {flag:>4} {r['chi2dof']:>10.3f} "
                      f"{r['chi2']:>9.2f} {r['ndof']:>5} {r['pvalue']:>6.3f} {tmc_s:>7}")

    summarise_selection(records, selection)

    def best_accepted_record(records, n, s, key="chi2dof"):
        """Accepted two-state record with smallest chi2/ndof for this (n, s); None if none."""
        acc = [r for r in records
               if r["nsquare"] == n and r["binsize"] == s and r["accepted"]
               and r["t_min_candidate"] is not None]
        return min(acc, key=lambda r: r[key]) if acc else None

    # ---- ONE-STATE GROUND-STATE EXTRACTION (uncorrelated) ----
    mass_records = []   # one row per (n, s)

    FIXED_WINDOW = (18, 39)

    FIXED_WINDOW_BY_ENSEMBLE = {"D251": (18, 39), "A650": (10, 21)}
    FIXED_WINDOW = FIXED_WINDOW_BY_ENSEMBLE[ensemble]

    for n in momentum_shells:
        for s in binsizes:
            rec = best_accepted_record(records, n, s)
            if rec is None:
                print(f"nsquare={n} binsize={s}: no accepted two-state, skip one-state")
                continue

            # t_min, t_max = rec["t_min_candidate"], rec["t_max"]
            t_min, t_max = FIXED_WINDOW 

            c2pt_per_nsquare = load_c2pt_per_nsquare(ensemble, momentum_shells, s)
            jkn = c2pt_per_nsquare[n]

            params_start = exclude_keys(
                estimate_c2pt_starting_values(c2pt_jkn=jkn, t_start=t_min, t_final=t_max),
                ["A1", "dE1"],
            )
            params_limit = {"A0": (0, None), "E0": (0, None)}

            fit_spec = {
                **base_spec,
                "binsize":              s,
                "correlation_type":     "uncorrelated",
                "nsquare":              n,
                "execute_resample_fit": True,        # needed for jackknife error on E0
                "model_id":             "one-state-exp",
                "t_start":              t_min,
                "t_final":              t_max,
                "params_start":         params_start,
                "params_limit":         params_limit,
            }

            one_state_fit = run_c2pt_fits(
                t_start              = fit_spec["t_start"],
                t_final              = fit_spec["t_final"],
                y                    = jkn,
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
            result   = fit_results(
                run_id          = run.run_id,
                fit_id          = fit_id,
                fit_hash        = fit_hash,
                fit_spec        = fit_spec,
                central         = one_state_fit.central,
                resample        = one_state_fit.resample,
                resample_method = fit_spec["resample_type"],
            )
            path_one = save_fit(result, run_dir, fit_id)
            print_fit_summary("one-state", one_state_fit, path=path_one)

            c      = one_state_fit.central
            E0     = float(c.params_est["E0"])
            E0_err = jackknife_err(one_state_fit, "E0")
            mass_records.append({
                "nsquare": n, "binsize": s, "t_min": t_min, "t_max": t_max,
                "E0": E0, "E0_err": E0_err,
                "chi2": float(c.chi2), "ndof": int(c.ndof),
                "chi2dof": float(c.chi2 / c.ndof) if c.ndof > 0 else float("inf"),
                "pvalue": float(c.pvalue),
                "src_t_zero": rec["t_zero"], "src_chi2dof": rec["chi2dof"],
            })

    # ---- BINSIZE DEPENDENCE OF THE NUCLEON MASS ----
    print("\n" + "=" * 84)
    print("NUCLEON MASS vs BINSIZE  (one-state, uncorrelated)")
    print("=" * 84)
    print(f" {'bin':>3} {'window':>9} {'E0':>12} {'E0_err':>10} {'chi2/ndof':>10} {'p':>6}  {'src t0':>6}  gvar")
    for r in sorted(mass_records, key=lambda r: (r["nsquare"], r["binsize"])):
        win = f"[{r['t_min']},{r['t_max']}]"
        print(f" {r['binsize']:>3} {win:>9} {r['E0']:12.6e} {r['E0_err']:10.3e} "
              f"{r['chi2dof']:10.3f} {r['pvalue']:6.3f}  {r['src_t_zero']:>6}  {gv.gvar(r['E0'], r['E0_err'])}")

    # mass_report = {
    #     "ensemble":        ensemble,
    #     "momentum_shells": momentum_shells,
    #     "binsizes":        binsizes,
    #     "selection_rule":  "lowest chi2/ndof accepted two-state per binsize",
    #     "fit":             {"model": "one-state-exp", "correlation": "uncorrelated"},
    #     "mass_records":    mass_records,
    # }
    # mass_path = run_dir / "nucleon_mass_vs_binsize.json"
    # with open(mass_path, "w") as f:
    #     json.dump(mass_report, f, indent=2)
    # print(f"\nwrote {mass_path}")
