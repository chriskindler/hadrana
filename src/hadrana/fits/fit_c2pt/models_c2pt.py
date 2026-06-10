# hadrana/fits/fit_c2pt/models.py
import numpy as np

PARAMS_ONE_STATE = {"A0", "E0"}
PARAMS_TWO_STATE = {"A0", "E0", "A1", "dE1"}

def model_c2pt_one_state(x, A0, E0):
    return A0 * np.exp(-E0 * x)

def model_c2pt_two_state(x, A0, E0, A1, dE1):
    return A0 * np.exp(-E0 * x) + A1 * np.exp(-(E0 + dE1) * x)

MODEL_REGISTRY = {
    "one-state-exp": (model_c2pt_one_state, PARAMS_ONE_STATE),
    "two-state-exp": (model_c2pt_two_state, PARAMS_TWO_STATE),
}

def get_model(model_id: str):
    """Return (model_callable, set of log-parametrized parameters)."""
    if model_id not in MODEL_REGISTRY:
        raise ValueError(
            f"Unknown model_id={model_id!r}. Valid: {list(MODEL_REGISTRY)}"
        )
    return MODEL_REGISTRY[model_id]

def estimate_c2pt_starting_values(c2pt_jkn, t_start, t_final):
    c2pt_avg = np.mean(c2pt_jkn, axis=0)
    t_eff    = (t_start + t_final) // 2
    E0_start = float(np.log(c2pt_avg[t_eff] / c2pt_avg[t_eff + 1]))
    A0_start = float(c2pt_avg[t_eff] * np.exp(E0_start * t_eff))
    return {
        "A0":  A0_start,
        "E0":  E0_start,
        "A1":  A0_start,
        "dE1": 0.25,
    }