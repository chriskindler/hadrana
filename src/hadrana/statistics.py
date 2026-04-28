import numpy as np

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
        Generate jackknife resampled ensembles from original correlator data (input_data) and corresponding reweighting factors (rwfs)
        Note:
        Jackknife ensembles: n_jkn  = n_cfg
        Bootstrap ensembles: n_bst != n_cfg
    """
    
    data, rwfs = bin_data(input_data, rwfs, bin_size)

    # axis=0 of any input data is always the number of configurations
    n_cfg = data.shape[0]
    n_rwf = rwfs.shape[0]

    data_cen = np.mean(data, axis=0) # overline{wX}
    rwfs_cen = np.mean(rwfs, axis=0) # overline{w}

    # jackknife resamples: (wO)^(J)_j and w^(J)_j
    data_jkn = (n_cfg * data_cen - data) / (n_cfg - 1)
    rwfs_jkn = (n_cfg * rwfs_cen - rwfs) / (n_cfg - 1)
    return (data_jkn / rwfs_jkn[:, np.newaxis, np.newaxis, np.newaxis])

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