import numpy as np 

def estimate_initial_timeslices(
    a_fm:         float,
    t_phys:       float,
    t_max:        int,
    t_start_phys: float = 0.13,   # physical earliest start [fm], ~lattice-2 on D251
    min_window:   int   = 6,      # two-state npar(4) + 2, so chigrad AICc is defined
) -> list[int]:
    """t_zero candidates in [t_start, t0_max], all leaving >= min_window
    points below t_max. Everything physical is fixed in fm; lattice indices
    are derived per ensemble via a_fm."""
    t_start = max(round(t_start_phys / a_fm), 1)
    t0_max  = round(t_phys / a_fm)

    # never start so late that [t_zero, t_max] is shorter than min_window
    t0_max = min(t0_max, t_max - min_window + 1)

    if t0_max < t_start:
        return []
    return list(range(t_start, t0_max + 1))

def estimate_maximum_timeslice(c2pt_jkn_avg, c2pt_jkn_err,
                               signal_to_noise_threshold,
                               boundary_margin=2):
    nt    = len(c2pt_jkn_avg)
    t_cut = nt // 2 - boundary_margin          # stay clear of the turn-up
    snr   = c2pt_jkn_avg / c2pt_jkn_err

    # first t where SNR drops below threshold AND stays below (sustained)
    below = snr[:t_cut] <= signal_to_noise_threshold
    tmax  = t_cut - 1
    for t in range(t_cut):
        if below[t] and below[t:].all():
            tmax = t - 1
            break
    return max(tmax, 0)

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
            t_min in [t_0, t_max] with A1 * exp(-dE1 * t) <= c2pt_err(t) / 4
            
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

if __name__ == "__main__":
    a = 0.064
    t_max = 44
    t_phys = 0.4
    t_start = 2
    hbarc = 0.1973

    timeslices = estimate_initial_timeslices(
        a_fm = a, 
        t_max = t_max,
        t_phys=t_phys,
        t_start=t_start,
        min_window=4
    )

    t_zero_max = hbarc/0.5
    t_phys_meas = t_zero_max / a
    print(t_phys_meas)
    print(f"t_zero_max = {hbarc / 0.5}")
    print(f"t_phys = {t_phys} fm")
    print(timeslices)
