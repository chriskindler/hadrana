import gvar as gv
import h5py
import pickle
import numpy as np

from typing import Literal, Optional

from chigrad.fit import (
    fit, FitResult,
    fit_constant_correlated, ConstantFitResult,
    compute_cov, compute_cov_inv,
)
from hadrana.fits.fit_c2pt.models_c2pt_ratio import (
    get_model,
    estimate_c2pt_ratio_starting_values,
)
from hadrana.fits.fit_c2pt.timeslice_criteria import (
    estimate_c2pt_ratio_minimum_timeslice,
    estimate_maximum_timeslice,
)

def perform_two_state_fit(
    t_min, t_max, y_res, cov_inv_full, start_params,
    tolerance=1e-4, strategy=1, ncall=10000,
):
    """Two-state Migrad fit on a window, using a sub-block of the full covariance."""
    model, _   = get_model("two-state")
    fit_range  = np.arange(t_min, t_max + 1)
    cov_inv    = cov_inv_full[np.ix_(fit_range, fit_range)]   # slice once
    y_window   = y_res[:, fit_range]

    return fit(
        t_data          = fit_range.astype(float),
        y_data          = y_window,
        correlated      = True,
        resample_method = "jackknife",
        resample_fit    = False,
        model           = model,
        start_params    = start_params,
        tolerance       = tolerance,
        strategy        = strategy,
        ncall           = ncall,
        cov_inv         = cov_inv,
    )

def perform_constant_fit(
    t_min, t_max, y_res, cov_inv_full, resample_fit,
):
    fit_range = np.arange(t_min, t_max + 1)
    cov_inv   = cov_inv_full[np.ix_(fit_range, fit_range)]    # slice once

    return fit_constant_correlated(
        t_data       = fit_range.astype(float),
        y_data       = y_res[:, fit_range],
        cov_inv      = cov_inv,
        resample_fit = resample_fit,
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

    log_w  = -chi2 / 2.0 + ndof
    w      = np.exp(log_w)
    p      = w / w.sum()

    mean = float(np.sum(cen * p))
    var  = float(np.sum((err**2 + cen**2) * p) - mean**2)
    error = float(np.sqrt(max(var, 0.0)))

    return mean, error, p

model_averages = []

def compute_exp_dEt_jkn(
    M_fit_jkn:  np.ndarray,
    nmax:    int,
    ns:     int,
    timeslices: np.ndarray
) -> np.ndarray:
    """
    E(p) = sqrt(p^2 + M^2)
    where p^2 = (2*pi/Ns)^2 * n^2
    dE(p) = sqrt(p^2 + M^2)
    """
    p_squared = (2.0 * np.pi / ns) ** 2 * nmax
    dE_jkn = np.sqrt(p_squared + M_fit_jkn ** 2) - M_fit_jkn
    return np.exp(dE_jkn[:, None] * timeslices[None, :])

# --------------------------------------------------------------------------------
# DRIVER
# --------------------------------------------------------------------------------

if __name__ == "__main__":
    ensemble           = "D251"
    nmax_values        = [1, 2, 3, 4, 5, 6, 8]
    initial_timeslices = [2, 3, 4, 5, 6, 7, 8, 9]
    snr_threshold      = 5

    ns = 64

    # ---- load data ----
    R_jkn_per_mom = {}
    path = f"/hdd/data/ensemble_data/{ensemble}/c2pt/{ensemble}_c2pt_jkn.h5"
    with h5py.File(path, "r") as file:
        for nmax in nmax_values:
            R_jkn_per_mom[nmax] = file[f"/c2pt_ratio/nsquare_{nmax}"][()]

    with open(
        f"/hdd/data/fit_results/{ensemble}/nucleon_energies/"
        f"{ensemble}_E0_jkn_mom00_tzero03.pkl", "rb"
    ) as f:
        M_fit_jkn = pickle.load(f)

    # ---- analyze one momentum ----
    for nmax in nmax_values:
        print(f"+++++++++++++++++++++ nmax = {nmax} ++++++++++++++++++++++++++")
        R_jkn       = R_jkn_per_mom[nmax]
        n_res, n_t  = R_jkn.shape
        timeslices  = np.arange(n_t)

        R_jkn       = R_jkn * compute_exp_dEt_jkn(M_fit_jkn, nmax, ns, timeslices)
        R_jkn_avg   = np.mean(R_jkn, axis=0)
        R_jkn_err   = np.std(R_jkn, axis=0, ddof=0) * np.sqrt(n_res - 1)

        t_max = estimate_maximum_timeslice(
            c2pt_jkn_avg              = R_jkn_avg,
            c2pt_jkn_err              = R_jkn_err,
            signal_to_noise_threshold = snr_threshold,
        )
        print(f"t_max = {t_max}")

        cov_full     = compute_cov(R_jkn, resample_method="jackknife")
        cov_inv_full = compute_cov_inv(cov_full)
        print(f"covariance shape: {cov_full.shape}")

        t0_list, tmin_list       = [], []
        B0_cen_list, B0_err_list = [], []
        chi2_list, ndof_list     = [], []

        for t0 in initial_timeslices:
            # two-state diagnostic fit
            start_params = estimate_c2pt_ratio_starting_values(
                c2pt_ratio_exp_jkn = R_jkn,
                initial_timeslice  = t0,
                maximum_timeslice  = t_max,
            )
            central_two_state, _ = perform_two_state_fit(
                t_min=t0, t_max=t_max, y_res=R_jkn,
                cov_inv_full=cov_inv_full,
                start_params=start_params,
            )
            v   = central_two_state.values
            B0_ts, B1_ts, dE1_ts = v["B0"], v["B1"], v["dE1"]

            # t_min from excited-state cutoff
            t_min = estimate_c2pt_ratio_minimum_timeslice(
                c2pt_ratio_jkn_err = R_jkn_err,
                initial_timeslice  = t0,
                maximum_timeslice  = t_max,
                B0_est             = B0_ts,
                B1_est             = B1_ts,
                dE1_est            = dE1_ts,
            )

            # reject degenerate two-state fits
            if t_min is None or dE1_ts > 1.0 or central_two_state.chi2_red > 2.5:
                print(f"  REJECTED t0={t0}: dE1={dE1_ts:.3f}, chi2/ndof={central_two_state.chi2_red:.2f}")
                continue

            # Production fit: closed-form weighted mean for B_0
            central_const, resample_const = perform_constant_fit(
                t_min=t_min, t_max=t_max, y_res=R_jkn,
                cov_inv_full=cov_inv_full,
                resample_fit=True,
            )

            B0_jkn  = np.array([r.values["B0"] for r in resample_const])
            B0_cen  = central_const.values["B0"]
            B0_err  = float(np.sqrt(n_res - 1) * np.std(B0_jkn, ddof=0))

            print(f"t0={t0:>2} tmin={t_min:>2}  "
                f"B0={gv.gvar(B0_cen, B0_err)}  "
                f"chi2/ndof={central_const.chi2_red:.2f}  ndof={central_const.ndof}")

            t0_list.append(t0)
            tmin_list.append(t_min)
            B0_cen_list.append(B0_cen)
            B0_err_list.append(B0_err)
            chi2_list.append(central_const.chi2)
            ndof_list.append(central_const.ndof)

        # ---- model average over t0 ----
        B0_avg, B0_avg_err, p = compute_model_average(
            np.asarray(chi2_list), np.asarray(ndof_list),
            np.asarray(B0_cen_list), np.asarray(B0_err_list),
        )

        model_averages.append(gv.gvar(B0_avg, B0_avg_err))

        print(f"\nMODEL AVERAGE for n^2 = {nmax}")
        print(f"{'t0':>3} {'tmin':>4} {'p_i':>8}   {'B0':>20}   {'chi2/ndof':>10}")
        for i, t0 in enumerate(t0_list):
            chi2_red = chi2_list[i] / ndof_list[i]
            print(f"{t0:>3d} {tmin_list[i]:>4d} {p[i]:>8.4f}   "
                f"{str(gv.gvar(B0_cen_list[i], B0_err_list[i])):>20}   "
                f"{chi2_red:>10.2f}")
        #print(f"\nB_0(n^2={nmax}) = {gv.gvar(B0_avg, B0_avg_err)}")

        for nmax, mavg in zip(nmax_values, model_averages):
            print(f"{nmax} & {mavg}")

