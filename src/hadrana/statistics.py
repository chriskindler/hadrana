import numpy as np

def maximum_binsize(ncfg: int):
    # magic number: 100 = number of bins
    return ncfg // 100

def generate_binsize_interval(maximum_binsize: int) -> list[int]:
    FINE_MAX   = 20    # dense integer block: where the rise of R(S) lives
    COARSE_MAX = 100   # sparse tail: flat region, adjacent points redundant
    COARSE_STEP = 10

    binsizes = []

    # fine region: every integer 1 .. min(20, s_max)
    for s in range(1, min(FINE_MAX, maximum_binsize) + 1):
        binsizes.append(s)

    # coarse tail: 30, 40, ... up to min(100, s_max)
    s = FINE_MAX
    while s + COARSE_STEP <= min(COARSE_MAX, maximum_binsize):
        s += COARSE_STEP
        binsizes.append(s)

    # always anchor the ceiling
    binsizes.append(maximum_binsize)

    return sorted(set(binsizes))

def bin_data(data: np.ndarray, rwfs: np.ndarray, bin_size: int):
    """
    input:
        data: c2pt/c3pt input data (n_cfg, n_src, n_mom, n_t)
        rwfs: reweighting factors  (n_cfg,)
        bin_size: bin size (no binning corresponds to bin_size = 1)
    output:
        sum (rwfs * data) / bin_size
        sum  rwfs / bin_size
    """
    n_cfg = data.shape[0]
    n_rwf = rwfs.shape[0]

    if n_cfg != n_rwf:
        raise RuntimeError(f"Number of configs ({n_cfg}) does not match number of reweighting factors ({n_rwf}).")

    # compute number of bins
    # trailing configs after flooring are discarded
    n_bin = n_cfg // bin_size
    
    data_bin = np.zeros((n_bin,) + data.shape[1:], dtype=np.float64)
    rwfs_bin = np.zeros((n_bin,)                 , dtype=np.float64)

    for i in range(n_bin):
    
        rwfs_data_sum = np.zeros(data.shape[1:], dtype=np.float64)
        rwfs_sum      = 0.0
        # range of j involves no (i + 1) * bin_size - 1 since python is already exclusive on the upper end:
        # range(a,b) gives the list(=interval) [a, a+1, ..., b-1]
        for j in range(bin_size * i, (i + 1) * bin_size):
            rwfs_data_sum += rwfs[j] * data[j]
            rwfs_sum      += rwfs[j]

        data_bin[i] = rwfs_data_sum / bin_size 
        rwfs_bin[i] = rwfs_sum / bin_size
    return data_bin, rwfs_bin 

def generate_jackknife_resamples(input_data: np.ndarray, rwfs: np.ndarray, bin_size: int):
    """
    Generate jackknife resampled ensembles from original correlator data (input_data)
    and corresponding reweighting factors (rwfs).

    Note:
        Jackknife ensembles: n_jkn  = n_cfg
        Bootstrap ensembles: n_bst != n_cfg
    """
    data, rwfs = bin_data(input_data, rwfs, bin_size)
    n_cfg = data.shape[0]

    data_cen = np.mean(data, axis=0)  # \overline{wX}
    rwfs_cen = np.mean(rwfs, axis=0)  # \overline{w}

    data_jkn = (n_cfg * data_cen - data) / (n_cfg - 1)
    rwfs_jkn = (n_cfg * rwfs_cen - rwfs) / (n_cfg - 1)

    # Reshape rwfs_jkn to (n_cfg, 1, 1, ...) matching the rank of data.
    rwfs_jkn = rwfs_jkn.reshape((n_cfg,) + (1,) * (data.ndim - 1))

    return data_jkn / rwfs_jkn

# TODO: Include option of binning for bootstrap construction
def generate_bootstrap_ensemble(data: np.ndarray, rwfs: np.ndarray, n_bst: int, bst_sample_index: int):
    n_cfg = data.shape[0]
    if rwfs.shape[0] != n_cfg:
        raise RuntimeError(f"Mismatch between #configs data ({n_cfg}) and rwfs ({rwfs.shape[0]}).")

    data_bst = np.zeros((n_bst, *data.shape[1:]), dtype=np.result_type(data, rwfs))

    for b in range(n_bst):
        bidx = bst_sample_index[b]
        y = data[bidx]
        w = rwfs[bidx].reshape(n_cfg, *([1] * (data.ndim - 1)))
        data_bst[b] = np.sum(w * y, axis=0) / np.sum(w)
    return data_bst

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
