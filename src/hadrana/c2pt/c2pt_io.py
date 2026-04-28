"""
data generator for c2pt data averaged over momenta and source ids

HDF5 file structure:

D251_c2pt_jkn.h5/
    nsquare_0/
        /fwd
            /source_set0
            /source_set1
        /fwd_bwd_avg
    nsquare_1/
        /fwd
            /source_set0
            /source_set1
        /fwd_bwd_avg
    nsquare_2/
        /fwd
            /source_set0
            /source_set1
        /fwd_bwd_avg
    nsquare_3/
        /fwd
            /source_set0
            /source_set1
        /fwd_bwd_avg
    nsquare_4/
        /fwd
            /source_set0
            /source_set1
        /fwd_bwd_avg
    nsquare_5/
        /fwd
            /source_set0
            /source_set1
        /fwd_bwd_avg
    nsquare_6/
        /fwd
            /source_set0
            /source_set1
        /fwd_bwd_avg
    nsquare_8/
        /fwd
            /source_set0
            /source_set1
        /fwd_bwd_avg

def get_source_set_match(ensemble: str):
    # save c2pt_fwd_jkn data from source set that must be combined with c3pt measurements
    pass

D251_c2pt_ratio_jkn.h5
    nsquare_1/
    nsquare_2/
    nsquare_3/
    nsquare_4/
    nsquare_5/
    nsquare_6/
    nsquare_8/
"""
