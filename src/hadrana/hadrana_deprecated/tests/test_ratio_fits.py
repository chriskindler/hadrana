import pickle
import numpy as np
import gvar as gv
import argparse

import iminuit
from iminuit import Minuit

from hadrana.loader import *
from hadrana.fits import *

if __name__ == "__main__":
    ensemble_id = "D251"
    source_sink_separations = [11, 14, 16, 19]
    fit_range = [1,1]

    ensemble_data_dict = load_bootstrap_data(ensemble_id, renormalised = True)

    # prepare xvals for combined ratio fit over all source sink separations
    # x is the set of all current insertion tau timeslices for all source-sink separations t
    x = np.array([[t, tau] for t in source_sink_separations for tau in range(0 + fit_range[0], t + 1 - fit_range[1])])

    print(ensemble_data_dict.keys())

    print(x)