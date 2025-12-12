import gvar as gv

"""
TODO
Must be transformed into a class including kappas, lattice spacing, source-sink separations etc.
"""

# ONLY THREE ENSEMBLES PER TRAJECTORY FOR NOW 
ENSEMBLE_IDS_SYM = ["D250", "X251", "X250"]
ENSEMBLE_IDS_TRM = ["D200", "D251", "N203"]
ENSEMBLE_IDS_MSC = ["D201", "N201", "N204"]

SOURCE_SINK_SEPARATIONS = {
    "SYM": {
        "X250": [11, 14, 16, 19],
        "X251": [11, 14, 16, 19],
        "D250": [11, 14, 16, 19],
    },
    "TRM": {
        "N203": [11, 14, 16, 19],
        "D251": [11, 14, 16, 19],
        "D200": [11, 14, 16, 19],
    },
    "MSC": {
        "N204": [11, 14, 16, 19],
        "N201": [11, 14, 16, 19],
        "D201": [11, 14, 16, 19],
    }
} 

# MESON MASSES IN ORDER OF ASCENDING PION MASSES PER TRAJECTORY 
HADRON_MASS_MAP = {
    "pion": {
        # SYM
        "D250": gv.gvar("0.06401(22)"),
        "X251": gv.gvar("0.08678(40)"),
        "X250": gv.gvar("0.11319(39)"),
        # TRM 
        "D200": gv.gvar("0.06540(33)"),
        "D251": gv.gvar("0.09203(16)"),
        "N203": gv.gvar("0.11245(30)"),
        # MSC
        "D201": gv.gvar("0.06472(42)"),
        "N201": gv.gvar("0.09268(31)"),
        "N204": gv.gvar("0.11423(33)"),
    },
    "kaon": {
        # SYM
        "D250": gv.gvar("0.06401(22)"),
        "X251": gv.gvar("0.08678(40)"),
        "X250": gv.gvar("0.11319(39)"),
        # TRM 
        "D200": gv.gvar("0.15652(15)"),
        "D251": gv.gvar("0.1503(21)"),
        "N203": gv.gvar("0.14399(24)"),
        # MSC
        "D201": gv.gvar("0.16302(18)"),
        "N201": gv.gvar("0.17040(22)"),
        "N204": gv.gvar("0.17734(29)"),
    },
    "nucleon": {
        # SYM
        "D250": gv.gvar("0.2886(28)"),
        "X251": gv.gvar("0.3182(85)"),
        "X250": gv.gvar("0.3594(51)"),
        # TRM 
        "D200": gv.gvar("0.3156(17)"),
        "D251": gv.gvar("0.3201(50)"),
        "N203": gv.gvar("0.3624(18)"),
        # MSC
        "D201": gv.gvar("0.3201(50)"),
        "N201": gv.gvar("0.3415(32)"),
        "N204": gv.gvar("0.3723(24)"),
    } 
}

def get_spatial_extent(id: str):
    spatial_extents = {
        "X250": 48, "X251": 48, "D250": 64, # SYM
        "N203": 48, "D251": 64, "D200": 64, # TRM
        "N204": 48, "N201": 48, "D201": 64  # MSC

    }
    spatial_extent = spatial_extents[id]
    return spatial_extent

def get_trajectory(id: str):
    if id in ["X250", "X251", "D250"]:
        trajectory = "SYM"
    elif id in ["N203", "D251", "D200"]:
        trajectory = "TRM"
    elif id in ["N204", "N201", "D201"]:
        trajectory = "MSC"
    return trajectory

def get_hadron_mass(id: str, hadron_mass_map: dict, hadron: str):
    hadron_mass = hadron_mass_map[hadron][id]
    return hadron_mass