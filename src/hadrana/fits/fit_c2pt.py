import argparse
import h5py
import json
import numpy as np

from datetime import datetime, timezone

from hadrana.loader import load_c2pt_jkn
from cora.fits.fit import (
    jackknife_covariance_matrix,
    inverse_covariance_matrix,
    CorrelatedChi2,
    minimise
)

"""
    Vary initial fit timeslice (1, 2, 3, 4, 5, 6)
    Take timeslice with best chi2 value?
    t_min_selection x {one_state, two_state}
"""

def build_c2pt_fits_filename(fit_params: dict) -> str:
    when = datetime.now(timezone.utc)
    time = when.strftime("%d-%m-%Y-%H-%M")

    if fit_params["resample_type"] == "jackknife":
        resample_type = "jkn"
    else:
        resample_type = "bst"

    format = "h5"

    ensemble      = fit_params["ensemble"]
    bin_size      = fit_params["bin_size"]
    stn_threshold = fit_params["c2pt"]["signal_to_noise_threshold"]
    tmin          = fit_params["c2pt"]["minimum_timeslice"]

    return f"{ensemble}-C2PT-FITS-{resample_type.upper()}-BINSIZE{bin_size}-STN{stn_threshold}-TSTART{tmin}-{time}.{format}"

def load_c2pt_jkn(ensemble: str):
    pass

def export_c2pt_fits(ensemble: str):
    pass

# --------------------------------------------------------------------------------
# MODELS
# --------------------------------------------------------------------------------

def model_c2pt_one_state(x, b0, E0):
    # b0 is log-amplitude; exp() keeps Migrad parameter O(1) and amplitude > 0
    b0 = np.exp(b0)
    return b0 * np.exp(-E0 * x)


def model_c2pt_two_state(x, b0, b1, E0, dE1):
    # b0, b1 are log-amplitudes; dE1 is log-gap
    # E1 = E0 + exp(dE1) guarantees E1 > E0 without needing Migrad limits
    b0, b1 = np.exp(b0), np.exp(b1)
    E1 = E0 + np.exp(dE1)
    return b0 * np.exp(-E0 * x) + b1 * np.exp(-E1 * x)

# --------------------------------------------------------------------------------
# PERFORM C2PT FIT 
# --------------------------------------------------------------------------------

def perform_c2pt_fit(c2pt: np.ndarray, nmax: int, fit_params: dict):
    n_res = c2pt.shape[0]

    SIGNAL_TO_NOISE_THRESHOLD: int = fit_params["c2pt"]["signal_to_noise_threshold"]
    MINIMUM_TIMESLICE:         int = fit_params["c2pt"]["minimum_timeslice"]

    c2pt_avg = np.mean(c2pt, axis=0)
    c2pt_err = np.std(c2pt, axis=0) * np.sqrt(n_res - 1)

    # ---- signal-to-noise based t_max (restricted to forward-propagating half) ----
    t_cut = c2pt.shape[1] // 2  # = 64 for Nt = 128
    StN_ratio = c2pt_avg / c2pt_err
    clean_timeslices = np.where(StN_ratio[:t_cut] > SIGNAL_TO_NOISE_THRESHOLD)[0]
    if len(clean_timeslices) == 0:
        raise RuntimeError(f"No clean timeslices below t={t_cut} at nmax={nmax}")
    t_max = int(clean_timeslices.max())

    # ---- initial fit range ----
    timeslices = np.arange(MINIMUM_TIMESLICE, t_max + 1)

    cov     = jackknife_covariance_matrix(c2pt[:, timeslices])
    cov_inv = inverse_covariance_matrix(cov, svd=False)

    # ---- starting values: physical -> log-parametrized ----
    start_phys = fit_params["c2pt"]["initial_parameters"][f"nmax-{nmax}"]
    start_log = {
        "A0":  np.log(start_phys["A0"]),  # log-amplitude
        "E0":  start_phys["E0"],
        "A1":  np.log(start_phys["A1"]),  # log-amplitude
        "dE1": np.log(start_phys["dE1"]), # log-gap: E1 = E0 + exp(dE1)
    }

    # ---- two-state fit (wide window, used to find t_min) ----
    two_state_chi2 = CorrelatedChi2(
        model   = model_c2pt_excited_state,
        x_data  = timeslices,
        y_data  = c2pt_avg[timeslices],
        cov_inv = cov_inv,
    )

    two_state_result = minimise(
        cost   = two_state_chi2,
        start  = start_log,
        limits = {"E0": (0, None)}, # exp() enforces dE1 > 0
    )

    log_params = {p: two_state_result.values[p] for p in start_log}

    A1_phys  = np.exp(log_params["A1"])
    dE1_phys = np.exp(log_params["dE1"])
    E1_phys  = log_params["E0"] + dE1_phys

    excited_state_contribution = A1_phys * np.exp(-E1_phys * timeslices)
    excited_state_exceeds      = excited_state_contribution <= (c2pt_err[timeslices] / 4.0)

    if np.any(excited_state_exceeds):
        t_min = int(timeslices[np.argmax(excited_state_exceeds)])
    else:
        print(
            f"[nmax={nmax}] WARNING: excited-state criterion never satisfied. "
            f"Falling back to t_min={MINIMUM_TIMESLICE}. "
            f"two_state valid={two_state_result.valid}, "
            f"E0_fit={log_params['E0']:.4f}, dE1_fit={dE1_phys:.4f}, E1_fit={E1_phys:.4f}"
        )
        t_min = MINIMUM_TIMESLICE

    # ---- one-state fit on refined window ----
    timeslices = np.arange(t_min, t_max + 1)
    cov     = jackknife_covariance_matrix(c2pt[:, timeslices])
    cov_inv = inverse_covariance_matrix(cov, svd=False)

    one_state_chi2 = CorrelatedChi2(
        model   = model_c2pt_ground_state,
        x_data  = timeslices,
        y_data  = c2pt_avg[timeslices],
        cov_inv = cov_inv,
    )

    one_state_result = minimise(
        cost   = one_state_chi2,
        start  = {"A0": log_params["A0"], "E0": log_params["E0"]},
        limits = {"E0": (0, None)},
    )

    log_central = {p: one_state_result.values[p] for p in ["A0", "E0"]}

    # ---- jackknife resample loop, same cov_inv ----
    A0_log_jkn = np.empty(n_res)
    E0_jkn     = np.empty(n_res)
    valid_jkn  = np.empty(n_res, dtype=bool)

    for j in range(n_res):
        cost_jkn = CorrelatedChi2(
            model   = model_c2pt_ground_state,
            x_data  = timeslices,
            y_data  = c2pt[j, timeslices],
            cov_inv = cov_inv,
        )

        result_jkn = minimise(
            cost   = cost_jkn,
            start  = log_central,
            limits = {"E0": (0, None)},
        )

        A0_log_jkn[j] = result_jkn.values["A0"]
        E0_jkn[j]     = result_jkn.values["E0"]
        valid_jkn[j]  = result_jkn.valid

    # apply exp() to each resample of the amplitude
    A0_jkn = np.exp(A0_log_jkn)

    # evaluate one-state fit function at resulting fit params per resample
    # shapes: A0_jkn (n_res,), E0_jkn (n_res,), timeslices (n_t,)
    # result: model_eval_jkn (n_res, n_t)
    model_eval_jkn = A0_jkn[:, None] * np.exp(-E0_jkn[:, None] * timeslices[None, :])

    n_density = 200
    timeslices_extended = np.linspace(t_min, t_max, n_density)
    model_eval_jkn_extended = A0_jkn[:, None] * np.exp(-E0_jkn[:, None] * timeslices_extended[None, :])

    ndof = one_state_result.ndof
    chi2 = one_state_result.chi2
    fmin = one_state_result.fmin

    return {
        "fit_range": (t_min, t_max),
        "fit_meta": {
            "valid":        one_state_result.valid,
            "converged":    one_state_result.converged,
            "accurate_cov": one_state_result.accurate_covariance,
            "fmin":         fmin,
        },
        "two_state_fit": {
            "valid": two_state_result.valid,
            "A0":    np.exp(log_params["A0"]),
            "A1":    A1_phys,
            "E0":    log_params["E0"],
            "dE1":   dE1_phys,
            "E1":    E1_phys,
        },
        "c2pt": c2pt,
        "timeslices": timeslices,
        "cov": cov,
        "cov_inv": cov_inv,
        "A0_jkn": A0_jkn, # physical amplitudes
        "E0_jkn": E0_jkn,
        "valid_jkn": valid_jkn,
        "timeslices_extended": timeslices_extended,
        "model_eval_jkn_extended": model_eval_jkn_extended,
        "ndof": ndof,
        "chi2": chi2,
        "chi2dof": chi2 / ndof,
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--filepath", help="json file path")
    args = parser.parse_args()

    with open(args.filepath, "r") as f:
        fit_params: dict = json.load(f)

    nmax_values = [0, 1, 2, 3, 4, 5, 6, 8]
    ensemble    = fit_params["ensemble"]
    source_sets = fit_params["source_sets"]
    BIN_SIZE    = fit_params["bin_size"]

    exceptionals = get_exceptional_indices(ensemble)
    rwfs = load_rwfs(ensemble)
    c2pt = load_c2pt(ensemble, source_sets)
    rwfs = np.delete(rwfs, exceptionals, axis=0)
    c2pt = np.delete(c2pt, exceptionals, axis=0)

    c2pt_per_nmax = compute_c2pt_per_nmax(c2pt, rwfs, nmax_values, BIN_SIZE)

    fit_result = {}

    for nmax in nmax_values:
        data = c2pt_per_nmax[nmax]
        c2pt_fit = perform_c2pt_fit(data, nmax, fit_params)

        fit_result[nmax] = c2pt_fit

    filename = build_c2pt_fits_filename(fit_params=fit_params, format="pickle")

    filepath = "/Users/ck/projects/data/c2pt/results"

    export_c2pt_fits(fit_results=fit_result, filepath=filepath, filename=filename, format="pickle")
