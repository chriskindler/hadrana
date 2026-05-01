_EXCEPTIONAL_INDICES = {
    "D251": [482], # small rwf
}

def get_exceptional_indices(ensemble: str):
    return _EXCEPTIONAL_INDICES[ensemble]
