import argparse
import gvar as gv
import h5py
import json
import numpy as np

from iminuit.util import describe
from scipy.stats import chi2 as chi2_dist
from typing import Callable, Literal, Optional

from cora.fit import fit, FitResult
from hadrana.fits.fit_c2pt.models_c2pt import (
    get_model,
    estimate_c2pt_starting_values
)
from hadrana.fits.fit_c2pt.timeslice_criteria import (
    estimate_minimum_timeslice,
    estimate_maximum_timeslice,
)

# --------------------------------------------------------------------------------
# C2PT FIT WRAPPER
# --------------------------------------------------------------------------------

def perform_c2pt_fit(
    t_min:           int,
    t_max:           int,
    y_res:           np.ndarray,                              # (n_res, n_time)
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

    Returns
    -------
    central_result   : FitResult
    resample_results : list[FitResult] or None
    """
    model, params = get_model(model_id)

    # Slice the fit window once; cora.fit treats axis 1 as ndata.
    fit_range = np.arange(t_min, t_max + 1)
    t_data    = fit_range.astype(float)
    y_window  = y_res[:, fit_range]                           # (n_res, n_t_fit)

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
    )


# --------------------------------------------------------------------------------
# MODEL AVERAGE (Jay-Neil AIC; Eqs. 52, 54, 55 of arXiv:2309.05774)
# --------------------------------------------------------------------------------

def compute_model_average(
    chi2_list:    np.ndarray,
    n_dof_list:   np.ndarray,
    obs_cen_list: np.ndarray,
    obs_err_list: np.ndarray,
) -> tuple[float, float, np.ndarray]:
    """
    Jay-Neil AIC-weighted model average.

        log w_i = -chi2_i/2 + n_dof_i           (Eq. 52)
        p_i     = w_i / sum_j w_j               (Eq. 54)
        mean    = sum_i O_i p_i                 (Eq. 55)
        err^2   = sum_i (sigma_i^2 + O_i^2) p_i - mean^2

    Returns (mean, error, p_array).
    """
    chi2  = np.asarray(chi2_list,    dtype=float)
    ndof  = np.asarray(n_dof_list,   dtype=float)
    cen   = np.asarray(obs_cen_list, dtype=float)
    err   = np.asarray(obs_err_list, dtype=float)

    log_w  = -chi2 / 2.0 + ndof
    # log_w -= log_w.max()              # log-sum-exp stability
    w      = np.exp(log_w)
    p      = w / w.sum()

    mean = float(np.sum(cen * p))
    var  = float(np.sum((err**2 + cen**2) * p) - mean**2)
    error = float(np.sqrt(max(var, 0.0)))

    return mean, error, p

# --------------------------------------------------------------------------------
# DRIVER
# --------------------------------------------------------------------------------

if __name__ == "__main__":
    ensemble = "D251"
    nmax_values = [0, 1, 2, 3, 4, 5, 6, 8]
    initial_timeslices = [3, 4, 5, 6, 7, 8, 9, 10]
    snr_threshold = 5

    def exclude_keys(d, keys):
        return {x: d[x] for x in d if x not in keys}

    c2pt_jkn_per_mom = {}
    path = f"/hdd/data/ensemble_data/{ensemble}/c2pt/{ensemble}_c2pt_jkn.h5"
    with h5py.File(path, "r") as file:
        for nmax in nmax_values:
            c2pt_jkn_per_mom[nmax] = file[f"/c2pt/nsquare_{nmax}/fwd_bwd_avg"][()]
            #print(c2pt_jkn_per_mom[nmax].shape)

    c2pt_jkn     = c2pt_jkn_per_mom[1]
    n_res        = c2pt_jkn.shape[0]
    c2pt_jkn_avg = np.mean(c2pt_jkn, axis=0)
    c2pt_jkn_err = np.std(c2pt_jkn, axis=0, ddof=0) * np.sqrt(n_res - 1)

    t_max = estimate_maximum_timeslice(
        c2pt_jkn_avg=c2pt_jkn_avg,
        c2pt_jkn_err=c2pt_jkn_err,
        signal_to_noise_threshold=snr_threshold,
    )

    # --- per-t0 collectors for the model average ---
    t0_list      = []
    tmin_list    = []
    ndata_list   = []
    A0_cen_list  = []
    A0_err_list  = []
    E0_cen_list  = []
    E0_err_list  = []
    chi2_list    = []
    ndof_list    = []
    chi2red_list = []

    for t0 in initial_timeslices:
        # ---- two-state diagnostic fit to determine t_min ----
        start_params = estimate_c2pt_starting_values(
            c2pt_jkn=c2pt_jkn,
            initial_timeslice=t0,
            maximum_timeslice=t_max,
        )

        central_two_state, _ = perform_c2pt_fit(
            t_min=t0, t_max=t_max, y_res=c2pt_jkn,
            resample_method="jackknife", model_id="two-state",
            start_params=start_params, correlated=True, resample_fit=False,
        )

        v2 = central_two_state.values
        A1_phys  = np.exp(v2["A1"])
        dE1_phys = np.exp(v2["dE1"])

        t_min = estimate_minimum_timeslice(
            c2pt_jkn_err=c2pt_jkn_err,
            initial_timeslice=t0, maximum_timeslice=t_max,
            A1_est=A1_phys, E0_est=v2["E0"], dE1_est=dE1_phys,
        )

        # ---- one-state production fit on the cleaned window ----
        start_params_one_state = exclude_keys(start_params, {"A1", "dE1"})

        central_one_state, resample_one_state = perform_c2pt_fit(
            t_min=t_min, t_max=t_max, y_res=c2pt_jkn,
            resample_method="jackknife", model_id="one-state",
            start_params=start_params_one_state,
            correlated=True, resample_fit=True,
        )

        window = np.arange(t_min, t_max + 1)

        ndata = len(window)

        # central values (back-transform A0)
        v = central_one_state.values
        A0_phys_central = float(np.exp(v["A0"]))
        E0_central      = float(v["E0"])

        # jackknife errors from per-resample fits
        A0_log_jkn  = np.array([r.values["A0"] for r in resample_one_state])
        E0_jkn      = np.array([r.values["E0"] for r in resample_one_state])
        A0_phys_jkn = np.exp(A0_log_jkn)

        A0_err_jkn = float(np.sqrt(n_res - 1) * np.std(A0_phys_jkn, ddof=0))
        E0_err_jkn = float(np.sqrt(n_res - 1) * np.std(E0_jkn,      ddof=0))

        # collect for model average
        t0_list.append(t0)
        tmin_list.append(t_min)
        ndata_list.append(ndata)
        A0_cen_list.append(A0_phys_central)
        A0_err_list.append(A0_err_jkn)
        E0_cen_list.append(E0_central)
        E0_err_list.append(E0_err_jkn)
        chi2_list.append(central_one_state.chi2)
        ndof_list.append(central_one_state.ndof)
        chi2red_list.append(central_one_state.chi2_red)

        print(f"\n{20*'='} one-state fit, t_initial={t0}, t_max={t_max} {20*'='}")
        print(f"  fit window      = ({t_min}, {t_max})")
        print(f"  A0 (phys, jkn)  = {gv.gvar(A0_phys_central, A0_err_jkn)}")
        print(f"  E0       (jkn)  = {gv.gvar(E0_central,      E0_err_jkn)}")
        print(f"  chi2/ndof         = {central_one_state.chi2_red:.2f}  (ndof={central_one_state.ndof})")

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

    # p_A0 == p_E0 because the weights depend only on chi2 and ndof, not on
    # the observable. Use either for diagnostics.
    p = p_A0

    print(f"\n{20*'='} MODEL AVERAGE (Jay-Neil AIC) {20*'='}")
    print(f"{'t0':>3} {'tmin':>4} {'ndata':>4} {'p_i':>8}   {'A0':>22}   {'E0':>20}   {'chi2/ndof':>10}")
    for i, t0 in enumerate(t0_list):
        # print(f"{i} & {t0} & {tmin_list[i]} & {ndata_list[i]} & {p[i]:.3f} & {gv.gvar(A0_cen_list[i], A0_err_list[i])} &  {gv.gvar(E0_cen_list[i], E0_err_list[i])} & {chi2red_list[i]:.2f}")
        # print(
        #     f"{t0:>3d} {tmin_list[i]:>4d} {p[i]:>8.4f} {ndata_list[i]} "
        #     f"{str(gv.gvar(A0_cen_list[i], A0_err_list[i])):>22}   "
        #     f"{str(gv.gvar(E0_cen_list[i], E0_err_list[i])):>20}   "
        #     f"{chi2red_list[i]:>10.3f}"
        # )

        print(f" {t0} & {tmin_list[i]} & {gv.gvar(A0_cen_list[i], A0_err_list[i])} &  {gv.gvar(E0_cen_list[i], E0_err_list[i])} & {p[i]:.3f} & {ndata_list[i]} & {ndof_list[i]} & {chi2red_list[i]:.2f}")

    print(f"\nFinal model-averaged results:")
    print(f"  A0 = {gv.gvar(A0_avg, A0_avg_err)}")
    print(f"  E0 = {gv.gvar(E0_avg, E0_avg_err)}")



