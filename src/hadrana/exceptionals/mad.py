import numpy as np

def identify_small_rwf_indices(rwfs: np.ndarray) -> np.ndarray:
    """
        input: rwfs - np.ndarray: dimension (n_cfg,)
        function:
            Flag small rwfs if absolute value of particular rwf is smaller
            than 5% of the absolute value of the average over all rwfs
        output:
            Corresponding indices of rwfs that are flagged

    """
    small_rwf_indices = np.argwhere(np.abs(rwfs)/np.mean(np.abs(rwfs)) < 0.05)
    return small_rwf_indices

def identify_exceptionals(obs_jkn: np.ndarray, sigma_threshold: int) -> np.ndarray:
    """
    Identify exceptional configurations using the MAD on jackknife resamples.

    check for absolute(R - med(R)) / sigma] > 6
    
    obs_jkn: shape (n_cfg, ...) — reweighted jackknife resamples
    sigma_threshold: number of normalised MADs to flag (default 5)
    
    returns: 1D boolean array of shape (n_cfg,), True = exceptional
    """
    N_res = obs_jkn.shape[0]

    # median over configurations (axis=0)
    median = np.median(obs_jkn, axis=0)  # shape: (...)
    
    # absolute deviations from median
    abs_dev = np.absolute(obs_jkn - median)  # shape: (n_cfg, ...)
    
    # MAD: median of absolute deviations over configurations
    mad = np.median(abs_dev, axis=0)  # shape: (...)
    
    # normalised MAD (Gaussian consistency factor)
    sigma_mad = 1.4826 * mad * np.sqrt(N_res)
    
    # avoid division by zero where MAD is exactly 0
    sigma_mad = np.where(sigma_mad == 0, np.inf, sigma_mad)
    
    # flag as exceptional if any element across (...) exceeds sigma_threshold
    exceeds = abs_dev / sigma_mad > sigma_threshold  # shape: (n_cfg, ...)
    # flatten all alpha into one axis: (n_cfg, ...) -> (n_cfg, n_kinematic)
    exceeds_flat = exceeds.reshape(N_res, -1)
    # flag configuration j if it exceeds threshold at ANY alpha
    exceptional = np.any(exceeds_flat, axis=1)
    
    return exceptional
