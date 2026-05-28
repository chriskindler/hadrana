import numpy as np 

def estimate_maximum_timeslice(
    c2pt_jkn_avg: np.ndarray,
    c2pt_jkn_err: np.ndarray,
    signal_to_noise_threshold: int
) -> int:
    """
        Determined upper bound of fit range
        t_max in [t_0, t_cut] with SN > signal_to_noise_threshold
    """

    nt = len(c2pt_jkn_avg)
    t_cut = nt // 2

    # signal-to-noise ratio
    snr = c2pt_jkn_avg / c2pt_jkn_err

    # signal-to-noise cutoff criterion
    valid_timeslices = np.where(snr[:t_cut] > signal_to_noise_threshold)[0]

    maximum_timeslice = int(valid_timeslices.max())
    return maximum_timeslice

def estimate_c2pt_minimum_timeslice(
    c2pt_jkn_err:      np.ndarray,
    initial_timeslice: int,
    maximum_timeslice: int,
    A1_est:            float,
    E0_est:            float,
    dE1_est:           float
) -> int | None:
    """
        DESCRIPTION
            Determine lower bound of fit range
            t_min in [t_0, t_max] with A0 * exp(-dE1 * t) <= c2pt_err(t) / 4
            
        INPUT
            c2pt_jkn_err: jackknife errors shape (n_t,)  
            A0_est: Estimated excited-state amplitude from two-state fit 
            E1_est: Estimated excited-state energy gap dE1 = E1 - E0 from two-state fit
            initial_timeslice: Starting timeslice from which we initially start the fits 
            maximum_timeslice: Maximum timeslice determined from two-state fit

        RETURNS
            Minimum timeslice for resulting fit range
    """

    window = np.arange(initial_timeslice, maximum_timeslice + 1)
    E1_est = E0_est + dE1_est

    excited_state = A1_est * np.exp(-E1_est * window)

    threshold = c2pt_jkn_err[window] / 4.0
    satisfied = excited_state <= threshold

    if not np.any(satisfied):
        print("Excited-state cutoff criterion never satisfied.")
        return None

    minimum_timeslice = int(window[np.argmax(satisfied)])
    return minimum_timeslice


def estimate_c2pt_ratio_minimum_timeslice(
    c2pt_ratio_jkn_err: np.ndarray,
    initial_timeslice: int,
    maximum_timeslice: int,
    B0_est: float,
    B1_est: float,
    dE1_est: float,
) -> int | None:
    """
        DESCRIPTION
            Determine lower bound of fit range for the two-state ratio model
                f_1(t) = B_0 * (1 + B_1 * exp(-dE1 * t)).
            t_min in [t_0, t_max] is the smallest timeslice at which the
            excited-state contribution drops below one quarter of the
            statistical error:
                |B_0 * B_1 * exp(-dE1 * t)| <= c2pt_ratio_err(t) / 4.

        INPUT
            c2pt_ratio_jkn_err: jackknife errors of the corrected ratio,
                shape (n_t,)
            B0_est:  Estimated ground-state amplitude from two-state fit
            B1_est:  Estimated excited-state amplitude from two-state fit
            dE1_est: Estimated excited-state gap dE1 = E1 - E0 from
                     two-state fit
            initial_timeslice: Starting timeslice
            maximum_timeslice: Maximum timeslice determined from SNR criterion

        RETURNS
            Minimum timeslice for resulting fit range, or None if the
            excited-state cutoff criterion is never satisfied within the
            window.
    """
    window = np.arange(initial_timeslice, maximum_timeslice + 1)
    excited_state = np.abs(B0_est * B1_est) * np.exp(-dE1_est * window)
    threshold = c2pt_ratio_jkn_err[window] / 4.0
    satisfied = excited_state <= threshold
    if not np.any(satisfied):
        print("Excited-state cutoff criterion never satisfied.")
        return None
    minimum_timeslice = int(window[np.argmax(satisfied)])
    return minimum_timeslice