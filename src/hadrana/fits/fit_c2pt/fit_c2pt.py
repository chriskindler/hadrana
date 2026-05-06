import argparse
import gvar as gv
import h5py
import json
import numpy as np

from scipy.stats import chi2 as chi2_dist

from cora.fit import (
    jackknife_covariance_matrix,
    inverse_covariance_matrix,
    CorrelatedChi2,
    minimise
)

from hadrana.fits.fit_c2pt.timeslice_criteria import (
    estimate_minimum_timeslice,
    estimate_maximum_timeslice,
    estimate_c2pt_starting_values
)

# --------------------------------------------------------------------------------
# MODELS
# --------------------------------------------------------------------------------

def model_c2pt_two_state(x, A0, A1, E0, dE1):
    # A0, A0 are log-amplitudes; dE1 is log-gap
    # E1 = E0 + exp(dE1) guarantees E1 > E0 without needing Migrad limits
    #A0, A1 = np.exp(A0), np.exp(A1)
    #E1 = E0 + np.exp(dE1)
    E1 = E0 + dE1
    return A0 * np.exp(-E0 * x) + A1 * np.exp(-E1 * x)

def model_c2pt_one_state(x, A0, E0):
    # A0 is log-amplitude; exp() keeps Migrad parameter O(1) and amplitude > 0
    # A0 = np.exp(A0)
    return A0 * np.exp(-E0 * x)

# --------------------------------------------------------------------------------
# C2PT FIT DRIVER
# --------------------------------------------------------------------------------

def perform_c2pt_two_state_fit(
    c2pt_jkn: np.ndarray,
    initial_timeslice: int,
    signal_to_noise_threshold: int
) -> dict:
    """ 
        DESCRIPTION
            Perform two-state c2pt fit with total 
            two-state fit -> t_max -> one-state fit -> t_min -> one-state fit using [t_min, t_max]

        INPUT
            c2pt_jkn_avg: Average of jackknife resamples shape (nt,) 
            initial_timeslice: Starting timeslice t_0 from which we start the fits 
            fit_config: Imported dictionary of fit configuration

        RETURNS
            Dictionary results with keys
            "nmax"
            "initial_timeslice"
            "maximum_timeslice"
            "estimated_params": {
                "A0_est":
                "E0_est":
                "A1_est":
                "dE1_est":
            }
    """

    n_jkn, _ = c2pt_jkn.shape
    c2pt_jkn_avg = np.mean(c2pt_jkn, axis=0)
    c2pt_jkn_err = np.std(c2pt_jkn, axis=0, ddof=0) * np.sqrt(n_jkn - 1)

    t_max = estimate_maximum_timeslice(
        c2pt_jkn_avg,
        c2pt_jkn_err,
        signal_to_noise_threshold
    )

    # FIT RELATED
    # TODO: Implement starting value strategy "data_driven" and "heuristic"
    starting_values = estimate_c2pt_starting_values(c2pt_jkn_avg, initial_timeslice, t_max) 
    fit_range = np.arange(initial_timeslice, t_max + 1)

    #cov     = jackknife_covariance_matrix(c2pt_jkn[:, fit_range])
    cov = np.diag(c2pt_jkn_err[fit_range] ** 2)
    cov_inv = inverse_covariance_matrix(cov, svd=False)

    # CENTRAL VALUE CHI2
    two_state_cost = CorrelatedChi2(
        model   = model_c2pt_two_state,
        x_data  = fit_range,
        y_data  = c2pt_jkn_avg[fit_range],
        cov_inv = cov_inv
    )

    # CENTRAL VALUE TWO-STATE FIT
    two_state_result = minimise(
        cost  = two_state_cost,
        start = starting_values,
        limits = {"E0": (0, None), "A0": (0, None), "A1": (0, None), "dE1": (0, None)},


    )

    estimated_params = {p: two_state_result.values[p] for p in starting_values}

    A1_est_log  = estimated_params["A1"]
    E0_est      = estimated_params["E0"]
    dE1_est_log = estimated_params["dE1"]

    A1_est  = np.exp(A1_est_log)
    dE1_est = np.exp(dE1_est_log)

    t_min = estimate_minimum_timeslice(
        c2pt_jkn_err,
        initial_timeslice,
        t_max,
        A1_est_log,
        E0_est,
        dE1_est_log
    )

    if t_min is None:
        t_min = initial_timeslice

    return {
        "t_min": t_min,
        "t_max": t_max,
        "two_state_result": two_state_result
    }

def perform_c2pt_one_state_fit(
    c2pt_jkn: np.ndarray,
    initial_timeslice: int,
    signal_to_noise_threshold: int
) -> dict:
    """ 
        One-state central and resample fits on [t_min, t_max]  -> return result dict 
    """

    n_jkn, _ = c2pt_jkn.shape
    c2pt_jkn_avg = np.mean(c2pt_jkn, axis=0)
    c2pt_jkn_err = np.std(c2pt_jkn, axis=0, ddof=0) * np.sqrt(n_jkn - 1)

    two_state_result = perform_c2pt_two_state_fit(
        c2pt_jkn,
        initial_timeslice,
        signal_to_noise_threshold
    )

    t_min = two_state_result["t_min"]
    t_max = two_state_result["t_max"]
    two_state_estimates = two_state_result["two_state_result"]

    fit_range = np.arange(t_min, t_max + 1)

    # cov = jackknife_covariance_matrix(c2pt_jkn[:, fit_range])
    cov = np.diag(c2pt_jkn_err[fit_range] ** 2)
    cov_inv = inverse_covariance_matrix(cov, svd=False)

    # CENTRAL VALUE FIT
    one_state_cost = CorrelatedChi2(
        model = model_c2pt_one_state,
        x_data = fit_range,
        y_data = c2pt_jkn_avg[fit_range],
        cov_inv = cov_inv
    )

    starting_values = estimate_c2pt_starting_values(c2pt_jkn_avg, t_min, t_max)
    one_state_start = {k: starting_values[k] for k in ("A0", "E0")}

    one_state_result = minimise(
        cost = one_state_cost,
        start = one_state_start
    )

    estimated_params_cen = {p: one_state_result.values[p] for p in ["A0", "E0"]}

    # RESAMPLE FITS
    A0_jkn = np.empty(n_jkn)
    E0_jkn = np.empty(n_jkn)
    
    for j in range(n_jkn):
        one_state_cost_jkn = CorrelatedChi2(
            model = model_c2pt_one_state,
            x_data = fit_range,
            y_data = c2pt_jkn[j, fit_range],
            cov_inv = cov_inv
        )

        one_state_result_jkn = minimise(
            cost = one_state_cost_jkn,
            start = estimated_params_cen
        )

        A0_jkn[j] = one_state_result_jkn.values["A0"]
        E0_jkn[j] = one_state_result_jkn.values["E0"]

    #A0_jkn = np.exp(A0_jkn)

    A0_cen = np.mean(A0_jkn, axis=0)
    A0_err = np.std(A0_jkn, axis=0, ddof=0) * np.sqrt(n_jkn - 1)

    E0_cen = np.mean(E0_jkn, axis=0)
    E0_err = np.std(E0_jkn, axis=0, ddof=0) * np.sqrt(n_jkn - 1)

    A0 = gv.gvar(A0_cen, A0_err)
    E0 = gv.gvar(E0_cen, E0_err)

    # After the one-state fit, before building the return dict:
    chi2_val = one_state_result.chi2
    ndof     = one_state_result.ndof
    chi2_red = chi2_val / ndof if ndof > 0 else float("inf")
    pvalue   = float(chi2_dist.sf(chi2_val, ndof)) if ndof > 0 else 0.0

    return {
        # window
        "t_min":           t_min,
        "t_max":           t_max,

        # one-state fit
        "A0":              A0,
        "E0":              E0,
        "A0_jkn":          A0_jkn,
        "E0_jkn":          E0_jkn,

        # one-state quality
        "chi2":     chi2_val,
        "ndof":     ndof,
        "chi2_red": chi2_red,
        "pvalue":   pvalue,
        "valid":           one_state_result.valid,
        "converged":       one_state_result.converged,
        "accurate_cov":    one_state_result.accurate_covariance,

        "cov":             cov,
        "cov_inv":         cov_inv,
        "fit_range":       fit_range,

        # two-state diagnostics
        "two_state_A0":    two_state_estimates.values["A0"],
        "two_state_A1":    two_state_estimates.values["A1"],
        "two_state_E0":    two_state_estimates.values["E0"],
        "two_state_dE1":   two_state_estimates.values["dE1"],
        "two_state_chi2":  two_state_estimates.chi2,
        "two_state_ndof":  two_state_estimates.ndof,
        "two_state_valid": two_state_estimates.valid,
    }


def print_c2pt_fit_result(result: dict, ensemble: str, nmax: int, initial_timeslice: int) -> None:
    """One-line summary plus a couple of diagnostic lines."""
    valid_str = "OK " if result["valid"] and result["converged"] else "BAD"

    # gvar formats as "0.4198(21)" — central(error) — quite readable
    E0 = result["E0"]
    A0 = result["A0"]

    # Two-state diagnostics
    ts_E0  = result["two_state_E0"]
    ts_dE1 = result["two_state_dE1"]
    ts_E1  = ts_E0 + ts_dE1
    ts_valid_str = "OK " if result["two_state_valid"] else "BAD"

    print(
        f"{ensemble} nmax={nmax:>2d} t_0={initial_timeslice:>2d}  "
        f"| [{result['t_min']:>2d}, {result['t_max']:>2d}] "
        f"({result['t_max'] - result['t_min'] + 1:>2d} pts)  "
        f"| E0 = {E0}  A0 = {A0}  "
        f"| chi2/dof = {result['chi2_red']:5.2f} "
        f"(p = {result['pvalue']:5.3f}, ndof = {result['ndof']:>2f})  "
        f"| {valid_str}  "
        f"| 2-state: E1 = {ts_E1:.4f}, dE1 = {ts_dE1:.4f} {ts_valid_str}"
    )

def print_nmax_header(ensemble: str, nmax: int, n_t_0_values: int) -> None:
    print()
    print("=" * 100)
    print(f"  {ensemble}  nmax = {nmax}   ({n_t_0_values} t_0 values)")
    print("=" * 100)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--filepath", help="json file path")
    args = parser.parse_args()

    with open(args.filepath, "r") as f:
        fit_config: dict = json.load(f)

    ensemble                  = fit_config["ensemble"]["id"]
    initial_timeslices        = fit_config["c2pt"]["initial_timeslices"]
    nmax_values               = fit_config["c2pt"]["nmax_values"]
    signal_to_noise_threshold = fit_config["c2pt"]["signal_to_noise_threshold"]

    c2pt_jkn_per_nmax = {}
    path = f"/hdd/data/ensemble_data/{ensemble}/c2pt/{ensemble}_c2pt_jkn.h5"
    with h5py.File(path, "r") as file:
        for nmax in nmax_values:
            c2pt_jkn_per_nmax[nmax] = file[f"/c2pt/nsquare_{nmax}/fwd_bwd_avg"][()]

    for nmax in nmax_values:
        c2pt_jkn = c2pt_jkn_per_nmax[nmax]
        print_nmax_header(ensemble, nmax, len(initial_timeslices))

        arr = c2pt_jkn_per_nmax[nmax]
        avg = arr.mean(axis=0)
        err = arr.std(axis=0, ddof=0) * np.sqrt(arr.shape[0] - 1)
        print(f"\nnmax={nmax}: shape={arr.shape}, n_jkn={arr.shape[0]}, n_t={arr.shape[1]}")

        for t0 in initial_timeslices:
            try:
                result = perform_c2pt_one_state_fit(
                    c2pt_jkn                  = c2pt_jkn,
                    initial_timeslice         = t0,
                    signal_to_noise_threshold = signal_to_noise_threshold,
                )
            except Exception as e:
                print(f"{ensemble} nmax={nmax:>2d} t_0={t0:>2d}  |  FAILED: {type(e).__name__}: {e}")
                continue

            print_c2pt_fit_result(result, ensemble, nmax, t0)