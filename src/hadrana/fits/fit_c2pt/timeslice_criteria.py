import numpy as np 

def estimate_maximum_timeslice(
    c2pt_jkn_avg: np.ndarray,
    c2pt_jkn_err: np.ndarray,
    signal_to_noise_threshold: int
) -> int:
    """
        Determined upper bound of fit range
        t_max in [t_0, nt - sink_distance] with SN > signal_to_noise_threshold
    """

    nt = len(c2pt_jkn_avg)
    t_cut = nt // 2

    # signal-to-noise ratio
    snr = c2pt_jkn_avg / c2pt_jkn_err

    # signal-to-noise cutoff criterion
    valid_timeslices = np.where(snr[:t_cut] > signal_to_noise_threshold)[0]

    maximum_timeslice = int(valid_timeslices.max())
    return maximum_timeslice

def estimate_minimum_timeslice(
    c2pt_jkn_err: np.ndarray,
    initial_timeslice: int,
    maximum_timeslice: int,
    A1_est: float,
    E0_est: float,
    dE1_est: float
) -> int | None:
    """
        DESCRIPTION
            Determine lower bound of fit range
            t_min in [t_0, t_max] with A0 * exp(-dE1 * t) <= c2pt_err(t) / 4
            
        INPUT
            c2pt_jkn_err: jackknife errors shape (n_t,)  
            A0_est: Estimated excited-state amplitude from two-state fit 
            E1_est: Estimated excited-state energy gap dE1 = E1 - E0 from two-state fit
            initial_timeslice: Starting timeslice t_0 from which we start the fits 
            maximum_timeslice: Maximum timeslice t_max determined from two-state fit

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

def estimate_c2pt_starting_values(
    c2pt_jkn_avg: np.ndarray,
    initial_timeslice: int,
    maximum_timeslice: int
) -> dict[str, float]:
    """
        DESCRIPTION
            Estimate jackknife data driven starting values for two-state fit

        INPUT
            c2pt_jkn_avg: Average over jackknife resamples ndim (nt,)
            t_0: Initial timeslice  

        RETURNS
            Dictionary of estimated starting vales for
                Amplitudes A0, A0
                Ground-state energy E0

    """
    # approximate timeslice at which effective energy plateau occurs
    # initial fit range
    t_eff = (initial_timeslice + maximum_timeslice) // 2

    # ground-state
    E0_start = float(np.log(c2pt_jkn_avg[t_eff] / c2pt_jkn_avg[t_eff + 1]))
    A0_start = float(c2pt_jkn_avg[t_eff] * np.exp(E0_start * t_eff))

    # excited-state 10% of ground-state
    A1_start = 0.1 * A0_start
    dE1_start = 0.4

    return {
        "A0": np.log(A0_start),
        "E0": E0_start,
        "A1": np.log(A1_start),
        "dE1": np.log(dE1_start)
    }

def estimate_c2pt_ratio_starting_values(
    c2pt_jkn_avg: np.ndarray,
    initial_timeslice: int,
    maximum_timeslice: int
) -> dict[str, float]:
    """
        DESCRIPTION
            Estimate jackknife data driven starting values for two-state fit

        INPUT
            c2pt_jkn_avg: Average over jackknife resamples ndim (nt,)
            t_0: Initial timeslice  

        RETURNS
            Dictionary of estimated starting vales for
                Amplitudes A0, A0
                Ground-state energy E0

    """
    # approximate timeslice at which effective energy plateau occurs
    # initial fit range
    t_eff = (initial_timeslice + maximum_timeslice) // 2

    # ground-state
    E0_start = float(np.log(c2pt_jkn_avg[t_eff] / c2pt_jkn_avg[t_eff + 1]))
    A0_start = float(c2pt_jkn_avg[t_eff] * np.exp(E0_start * t_eff))

    # excited-state 10% of ground-state
    A1_start = 0.1 * A0_start
    dE1_start = 0.4

    return {
        "A0": np.log(A0_start),
        "E0": E0_start,
        "A1": np.log(A1_start),
        "dE1": np.log(dE1_start)
    }
