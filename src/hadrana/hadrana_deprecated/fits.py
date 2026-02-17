# DEFAULT
import numpy as np
import iminuit

def chi2_distribution(x, y, e, model: callable, *params,):
    y_fit = model(x, *params)
    delta = y - y_fit
    chi2 = (delta / e) ** 2
    return chi2


def generic_least_square(
    x: np.ndarray,
    y_cen: np.ndarray,
    y_err: np.ndarray,
    params: dict,
    priors: dict,
    model_id: str,
    model: callable
):
    """
    chi2 = chi2_data + chi2_prior
    chi2_data = sum_i (y_i - f(x_i, *params)) ** 2 / y_err_i ** 2
    chi2_prior = TBA

    params = {
        {"param_id": initial_values =  list[float]}
    }
    """

    x     = np.asarray(x,     dtype=float)
    y_cen = np.asarray(y_cen, dtype=float)
    y_err = np.asarray(y_err, dtype=float)

    param_ids = list(params[model_id].keys())
    param_val = list(params[model_id].values())

    chi2_data = ( (y_cen - model(x, *param_val) )  / y_err ) ** 2

    if priors is None:
        chi2 = chi2_data
    # TODO: add chi2 = chi2_data + chi2_prior

    return chi2

# CUSTOM
def fit(
    x: np.ndarray,  
    y: np.ndarray,
    params: dict,
    priors: dict,
    model: callable,
    correlated: bool = False
):
    pass
    


if __name__ == "__main__":

    import random
    import numpy as np

    params = {
        "B3": 
        {"S": 1.242, "b1": 0.242, "dM1": 0.232},
    }

    model_id = "B3"
  
    param_ids = list(params[model_id].keys())
    param_val = list(params[model_id].values())

    print(param_ids)
    print(param_val)





