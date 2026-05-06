# hadrana/fits/fit_c2pt/models.py
import numpy as np

# Convention: amplitudes A0, A1 are log-parametrized internally (Migrad sees logs);
# the gap dE1 is also log-parametrized (so E1 > E0 is guaranteed).
# E0 is in physical units.
# Conversion to physical happens in results.py.

PARAMS_ONE_STATE = {"A0", "E0"}
PARAMS_TWO_STATE = {"A0", "E0", "A1", "dE1"}

def model_c2pt_one_state(x, A0, E0):
    """One-state c2pt model. A0 is log-amplitude, E0 physical."""
    return np.exp(A0) * np.exp(-E0 * x)

def model_c2pt_two_state(x, A0, E0, A1, dE1):
    return np.exp(A0) * np.exp(-E0 * x) + np.exp(A1) * np.exp(-(E0 + np.exp(dE1)) * x)

MODEL_REGISTRY = {
    "one-state": (model_c2pt_one_state, PARAMS_ONE_STATE),
    "two-state": (model_c2pt_two_state, PARAMS_TWO_STATE),
}

def get_model(model_id: str):
    """Return (model_callable, set of log-parametrized parameters)."""
    if model_id not in MODEL_REGISTRY:
        raise ValueError(
            f"Unknown model_id={model_id!r}. Valid: {list(MODEL_REGISTRY)}"
        )
    return MODEL_REGISTRY[model_id]


def estimate_c2pt_starting_values(c2pt_jkn, initial_timeslice, maximum_timeslice):
    c2pt_avg = np.mean(c2pt_jkn, axis=0)
    t_eff = (initial_timeslice + maximum_timeslice) // 2
    E0_start  = float(np.log(c2pt_avg[t_eff] / c2pt_avg[t_eff + 1]))
    A0_phys   = float(c2pt_avg[t_eff] * np.exp(E0_start * t_eff))
    return {
        "A0":  np.log(A0_phys),           # log-amplitude, ~ -30
        "E0":  E0_start,                  # physical, ~ 0.34
        "A1":  np.log(0.1 * A0_phys),     # log-amplitude, ~ -32
        "dE1": np.log(0.4),               # log-gap, ~ -0.92
    }