SCALAR_KEYS = (
    # identification
    "fit_id",
    "fit_hash",
    "timestamp",
    "ensemble",
    "model_id",
    "resample_method",
    "correlation_type",
    "bin_size",
    # counts
    "ndata",
    "npar",
    "ndof",
    "nres",
    "nvalid",
    # quality
    "chi2",
    "chi2dof",
    "pvalue",
    "aic",
    "aicc",
    # diagnostic
    "valid",
    "converged",
    "accurate_cov",
    "fmin_edm",
    # meta
    "tolerance",
    "strategy",
    "ncall",
    # array bounds
    "fit_range_min",
    "fit_range_max",
)

ARRAY_KEYS = (
    "fit_range",
    "fit_range_ext",
    "residuals",
    "y_fit",
    "y_fit_cen_ext",
    "y_fit_err_ext",
    "valid_res",
)

DICT_KEYS = (
    "params_start",
    "params_limit",
    "params_cen",
    "params_err",
    "params_err_hesse",
    "params_res",
)