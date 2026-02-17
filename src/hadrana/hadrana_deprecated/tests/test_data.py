import numpy as np

from argparse import Namespace

import iminuit
from iminuit.cost import LeastSquares
from iminuit import Minuit, describe

from hadrana.loader import load_bootstrap_data
from hadrana.fits import generic_least_square

if __name__ == "__main__":

    ensemble_id = "D251"
    source_sink_separations = [11, 14, 16, 19]
    qmax = 0

    ensemble_data = load_bootstrap_data(ensemble_id, renormalised=True)

    test_key = f"{ensemble_id}-RATIO-REN-MOMENTUM-PERP-QMAX{qmax}-11"

    Nbst = 100

    y_cen = ensemble_data[ensemble_id][test_key][-1]
    y_bst = np.zeros( (Nbst, len(y_cen)))

    for b in range(Nbst):
        y_bst[b] = ensemble_data[ensemble_id][test_key][b]

    y_err = np.std(y_bst, axis=0, ddof=1)

    print(y_cen)
    print(y_err)


