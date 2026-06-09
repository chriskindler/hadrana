# hadrana/fits/fit_c2pt/models_c2pt_ratio.py
import numpy as np

# Convention: amplitudes A0, A1 are log-parametrized internally (Migrad sees logs);
# the gap dE1 is also log-parametrized (so E1 > E0 is guaranteed).
# E0 is in physical units.
# Conversion to physical happens in results.py.

PARAMS_CONSTANT = {"B0"}
PARAMS_TWO_STATE = {"B0", "B1", "dE1"}
PARAMS_TWO_STATE_MASS_GAP = {"B0", "B1", "dE1", "B10", "dM1"}

def weighted_mean_correlated(y: np.ndarray, cov_inv: np.ndarray):
    """
    Correlated weighted mean: minimizes (y - Z*1)^T C^-1 (y - Z*1).
    Returns Z, chi2, ndof.
    """
    ones = np.ones_like(y, dtype=float)
    denom = ones @ cov_inv @ ones
    Z     = (ones @ cov_inv @ y) / denom
    resid = y - Z
    chi2  = resid @ cov_inv @ resid
    ndof  = len(y) - 1
    return Z, chi2, ndof

def model_c2pt_ratio_constant(x, B0):
    return B0 * np.ones_like(x, dtype=float)

def model_c2pt_ratio_two_state(x, B0, B1, dE1):
    return B0 * (1 +  B1 * np.exp(-dE1 * x))

def model_c2pt_ratio_two_state_mass_gap(x, B0, B1, dE1, B10, dM1):
    return B0 * (1 +  B1 * np.exp(-dE1 * x) - B10 * np.exp(-dM1 * x) )

MODEL_REGISTRY = {
    "constant": (model_c2pt_ratio_constant, PARAMS_CONSTANT),
    "two-state": (model_c2pt_ratio_two_state, PARAMS_TWO_STATE),
    "two-state-mass-gap": (model_c2pt_ratio_two_state_mass_gap, PARAMS_TWO_STATE_MASS_GAP),
}

def get_model(model_id: str):
    """Return (model_callable, set of log-parametrized parameters)."""
    if model_id not in MODEL_REGISTRY:
        raise ValueError(
            f"Unknown model_id={model_id!r}. Valid: {list(MODEL_REGISTRY)}"
        )
    return MODEL_REGISTRY[model_id]

def estimate_c2pt_ratio_starting_values(
    c2pt_ratio_exp_jkn: np.ndarray,
    initial_timeslice: int,
    maximum_timeslice: int,
) -> dict[str, float]:
    """
    Starting values for the two-state ratio fit
        f_1(p, t) = B_0 * (1 + B_1 * exp(-dE1 * t)).
    """
    avg = np.mean(c2pt_ratio_exp_jkn, axis=0)

    B0_start  = float(avg[maximum_timeslice])
    B1_start  = 0.1
    dE1_start = 0.5

    return {
        "B0":  B0_start,
        "B1":  B1_start,
        "dE1": dE1_start,
    }