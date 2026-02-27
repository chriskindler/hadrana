import h5py
import numpy as np
import pickle
from collections import defaultdict

from cora.ensembles.ensemble_helpers import * 
from cora.reweighting import get_rwfs

# H106
# replica order
# replicas = ['r003', 'r001', 'r002','r004']
# source_set = ['source_set2', 'source_set4']

# source_set2:
# LHP_SEQ = { "10": [0,1], "12": [0,1], "14": [0,1]}
# source_set4:
# REG_SEQ = { "10": [0,1], "12": [0,1], "14": [0,1]}

# 1) Construct replica-concatinated forward-backward averaged c2pt
# 2) Construct replica-concatinated isovector axial and pseudoscalar c3pt
# 3) Construct ratios per source-sink distance and measurement
# 4) Bootstrap resmapling of ratios
# 5) Pickle and export bootstrap ratios

# replica order:
# "H106": [("r003",-1),("r001",-1),("r002",1),("r004",1)] 

if __name__ == "__main__":
    ensemble = "H106"
    replicas = ['r003', 'r001', 'r002', 'r004']
    
    print("Hello World")
