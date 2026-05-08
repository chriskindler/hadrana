import argparse
import gvar as gv
import h5py
import json
import pickle
import numpy as np

from iminuit.util import describe
from scipy.stats import chi2 as chi2_dist
from typing import Callable, Literal, Optional

from cora.fit import fit, FitResult, compute_cov, compute_cov_inv
from hadrana.fits.fit_c2pt.models_c2pt import (
    get_model,
    estimate_c2pt_starting_values,
)
from hadrana.fits.fit_c2pt.timeslice_criteria import (
    estimate_c2pt_minimum_timeslice,
    estimate_maximum_timeslice,
)

# --------------------------------------------------------------------------------
# C2PT FIT WRAPPER (with precomputed covariance)
# --------------------------------------------------------------------------------


def perform_c2pt_fit(
    t_min:           int,
    t_max:           int,
    y_res:           np.ndarray,                              # (n_res, n_t)
    cov_full:        np.ndarray,                              # (n_t, n_t) precomputed
    resample_method: Literal["bootstrap", "jackknife"],
    model_id:        Literal["one-state", "two-state"],
    start_params:    dict[str, float],
    correlated:      bool,
    resample_fit:    bool,
    limits:          Optional[dict[str, tuple[float | None, float | None]]] = None,
    tolerance:       float = 1e-4,
    strategy:        int   = 1,
    ncall:           int   = 10000,
) -> tuple[FitResult, Optional[list[FitResult]]]:
    """
    Thin wrapper around cora.fit.fit for C2PT one-/two-state fits.
    Uses a precomputed full covariance matrix and slices the relevant
    sub-block for the [t_min, t_max] fit window.
    """
    model, _ = get_model(model_id)

    fit_range = np.arange(t_min, t_max + 1)
    t_data    = fit_range.astype(float)
    y_window  = y_res[:, fit_range]                           # (n_res, n_t_fit)

    # Slice the sub-block of the full covariance and invert.
    if correlated:
        cov_window = cov_full[np.ix_(fit_range, fit_range)]
        cov_inv    = compute_cov_inv(cov_window)
        sdev_inv   = None
    else:
        cov_inv    = None
        # uncorrelated: use the diagonal of the sliced covariance for sdev
        sdev       = np.sqrt(np.diag(cov_full[np.ix_(fit_range, fit_range)]))
        sdev_inv   = 1.0 / sdev

    return fit(
        t_data          = t_data,
        y_data          = y_window,
        correlated      = correlated,
        resample_method = resample_method,
        resample_fit    = resample_fit,
        model           = model,
        start_params    = start_params,
        limits          = limits,
        tolerance       = tolerance,
        strategy        = strategy,
        ncall           = ncall,
        cov_inv         = cov_inv,
        sdev_inv        = sdev_inv,
    )

# --------------------------------------------------------------------------------
# MODEL AVERAGE (arXiv:2309.05774)
# --------------------------------------------------------------------------------

def compute_model_average(
    chi2_list:    np.ndarray,
    n_dof_list:   np.ndarray,
    obs_cen_list: np.ndarray,
    obs_err_list: np.ndarray,
) -> tuple[float, float, np.ndarray]:
    """
        AIC-weighted model average.
        log w_i = -chi2_i/2 + n_dof_i
        p_i     = w_i / sum_j w_j
        mean    = sum_i O_i p_i
        error^2 = sum_i (sigma_i^2 + O_i^2) p_i - mean^2

    Returns (mean, error, p_array).
    """
    chi2  = np.asarray(chi2_list,    dtype=float)
    ndof  = np.asarray(n_dof_list,   dtype=float)
    cen   = np.asarray(obs_cen_list, dtype=float)
    err   = np.asarray(obs_err_list, dtype=float)

    log_w = -chi2 / 2.0 + ndof
    w     = np.exp(log_w)
    p     = w / w.sum()

    mean  = float(np.sum(cen * p))
    var   = float(np.sum((err**2 + cen**2) * p) - mean**2)
    error = float(np.sqrt(max(var, 0.0)))

    return mean, error, p

# --------------------------------------------------------------------------------
# DRIVER
# --------------------------------------------------------------------------------

if __name__ == "__main__":
    ensemble    = "D251"
    nmax_values = [0, 1, 2, 3, 4, 5, 6, 8]
    initial_timeslices = [3, 4, 5, 6, 7]
    snr_threshold      = 5

    def exclude_keys(d, keys):
        return {x: d[x] for x in d if x not in keys}

    # ---- load data ----
    c2pt_jkn_per_mom = {}
    path = f"/hdd/data/ensemble_data/{ensemble}/c2pt/{ensemble}_c2pt_jkn.h5"
    with h5py.File(path, "r") as file:
        for nmax in nmax_values:
            c2pt_jkn_per_mom[nmax] = file[f"/c2pt/nsquare_{nmax}/fwd_bwd_avg"][()]

    for nmax in nmax_values:
        #nmax         = 0
        c2pt_jkn     = c2pt_jkn_per_mom[nmax]
        n_res        = c2pt_jkn.shape[0]
        c2pt_jkn_avg = np.mean(c2pt_jkn, axis=0)
        c2pt_jkn_err = np.std(c2pt_jkn, axis=0, ddof=0) * np.sqrt(n_res - 1)

        t_max = estimate_maximum_timeslice(
            c2pt_jkn_avg              = c2pt_jkn_avg,
            c2pt_jkn_err              = c2pt_jkn_err,
            signal_to_noise_threshold = snr_threshold,
        )

        cov_full = compute_cov(c2pt_jkn, resample_method="jackknife")
        print(f"covariance shape: {cov_full.shape}")

        t0_list      = []
        tmin_list    = []
        ndata_list   = []
        A0_cen_list  = []
        A0_err_list  = []
        E0_jkn_list  = []
        E0_cen_list  = []
        E0_err_list  = []
        chi2_list    = []
        ndof_list    = []
        chi2red_list = []

        for t0 in initial_timeslices:
            start_params = estimate_c2pt_starting_values(
                c2pt_jkn          = c2pt_jkn,
                initial_timeslice = t0,
                maximum_timeslice = t_max,
            )

            # TWO-STATE FIT
            central_two_state, _ = perform_c2pt_fit(
                t_min=t0, t_max=t_max, y_res=c2pt_jkn,
                cov_full=cov_full,
                resample_method="jackknife", model_id="two-state",
                start_params=start_params, correlated=True, resample_fit=False,
            )

            v2       = central_two_state.values
            A1_phys  = np.exp(v2["A1"])
            dE1_phys = np.exp(v2["dE1"])

            t_min = estimate_c2pt_minimum_timeslice(
                c2pt_jkn_err      = c2pt_jkn_err,
                initial_timeslice = t0,
                maximum_timeslice = t_max,
                A1_est            = A1_phys,
                E0_est            = v2["E0"],
                dE1_est           = dE1_phys,
            )

            # ---- one-state production fit ----
            start_params_one_state = exclude_keys(start_params, {"A1", "dE1"})

            central_one_state, resample_one_state = perform_c2pt_fit(
                t_min=t_min, t_max=t_max, y_res=c2pt_jkn,
                cov_full=cov_full,
                resample_method="jackknife", model_id="one-state",
                start_params=start_params_one_state,
                correlated=True, resample_fit=True,
            )

            window = np.arange(t_min, t_max + 1)
            ndata  = len(window)

            # central values (back-transform A0)
            v               = central_one_state.values
            A0_phys_central = float(np.exp(v["A0"]))
            E0_central      = float(v["E0"])

            # jackknife errors from per-resample fits
            A0_log_jkn  = np.array([r.values["A0"] for r in resample_one_state])
            E0_jkn      = np.array([r.values["E0"] for r in resample_one_state])
            A0_phys_jkn = np.exp(A0_log_jkn)

            A0_err_jkn = float(np.sqrt(n_res - 1) * np.std(A0_phys_jkn, ddof=0))
            E0_err_jkn = float(np.sqrt(n_res - 1) * np.std(E0_jkn,      ddof=0))

            # collect
            t0_list.append(t0)
            tmin_list.append(t_min)
            ndata_list.append(ndata)
            A0_cen_list.append(A0_phys_central)
            A0_err_list.append(A0_err_jkn)
            E0_jkn_list.append(E0_jkn)
            E0_cen_list.append(E0_central)
            E0_err_list.append(E0_err_jkn)
            chi2_list.append(central_one_state.chi2)
            ndof_list.append(central_one_state.ndof)
            chi2red_list.append(central_one_state.chi2_red)

            print(f"\n{20*'='} one-state fit, t_initial={t0}, t_max={t_max} {20*'='}")
            print(f"  fit window      = ({t_min}, {t_max})")
            print(f"  A0 (phys, jkn)  = {gv.gvar(A0_phys_central, A0_err_jkn)}")
            print(f"  E0       (jkn)  = {gv.gvar(E0_central,      E0_err_jkn)}")
            print(f"  chi2/ndof       = {central_one_state.chi2_red:.2f}  (ndof={central_one_state.ndof})")

        # ----------------------------------------------------------------
        # MODEL AVERAGE OVER t0
        # ----------------------------------------------------------------
        chi2_arr = np.asarray(chi2_list, dtype=float)
        ndof_arr = np.asarray(ndof_list, dtype=float)

        A0_avg, A0_avg_err, p_A0 = compute_model_average(
            chi2_arr, ndof_arr, np.asarray(A0_cen_list), np.asarray(A0_err_list),
        )
        E0_avg, E0_avg_err, p_E0 = compute_model_average(
            chi2_arr, ndof_arr, np.asarray(E0_cen_list), np.asarray(E0_err_list),
        )

        p = p_A0

        print(f"\n{20*'='} MODEL AVERAGE (AIC) {20*'='}")
        for i, t0 in enumerate(t0_list):
            print(
                f" {t0} & {tmin_list[i]} & "
                f"{gv.gvar(A0_cen_list[i], A0_err_list[i])} &  "
                f"{gv.gvar(E0_cen_list[i], E0_err_list[i])} & "
                f"{p[i]:.3f} & {chi2_list[i]:.2f} & {ndof_list[i]} & "
                f"{chi2red_list[i]:.2f}"
            )

        print(f"\n Model-averaged results:")
        print(f"  A0 = {gv.gvar(A0_avg, A0_avg_err)}")
        print(f"  E0 = {gv.gvar(E0_avg, E0_avg_err)}")

        base_path = f"/hdd/data/fit_results/{ensemble}/nucleon_energies"
        for i, t0 in enumerate(t0_list):
            file_name = f"{ensemble}_E0_jkn_mom0{nmax}_tzero0{t0}.pkl"
            E0_data   = E0_jkn_list[i]

            with open(f"{base_path}/{file_name}", "wb") as f:
                pickle.dump(E0_data, f, protocol=pickle.HIGHEST_PROTOCOL)