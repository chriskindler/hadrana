import gvar as gv
import h5py
import numpy as np
import pickle
import os

from cora.ensembles import EnsembleSpecs
from cora.exceptionals.mad import *
from cora.logging import log
from cora.momenta import * 
from cora.statistics import generate_jackknife_resamples

# TODO: Compare lattice vs. continuum energy
def compute_lattice_energy(ensemble: str, nmax: int, hadron: str):
    pass

def compute_continuum_energy(ensemble: str, nmax: int, hadron: str) -> np.ndarray:
    M           = get_hadron_mass(ensemble, "NUCLEON")
    Ns          = get_spatial_lattice_points(ensemble)
    n_vec_list  = get_selected_momenta(nmax)
    # q_j = (2pi / Ns) * n_j for each spatial direction
    # q^2 = sum_j q_j^2 is the same for all momenta in the shell
    n_vec            = np.array(n_vec_list) # lattice unit vectors, shape (N_mom, 3)
    q_vec            = (2 * np.pi / Ns) * n_vec # lattice momenta, shape (N_mom, 3)
    q_vec_squared    = np.sum(q_vec ** 2, axis=1) # square of lattice mommenta, shape (N_mom,)
    continuum_energy = np.sqrt(q_vec_squared + M ** 2)
    return continuum_energy

def compute_jackknife_c2pt_ratio(
    ensemble:   str,
    source_set: str,
    c2pt:       np.ndarray,
    rwfs:       np.ndarray,
    nmax:       int,
    bin_size:   int) -> np.ndarray:

    # from 93 available momentum configurations nvec = (nx, ny, nz)
    # we select all nvec that satisfy nx^2 + ny^2 + nz^2 = nmax
    # we always consider unpolarised two-point functions
    mom_indices: list[tuple(int, int, int)] = get_selected_momenta_indices(nmax, pol=None)

    c2pt_q = c2pt[:, :, mom_indices, :]
    c2pt_0 = c2pt[:, :, 0, :]
    c2pt_0 = c2pt_0[:, :, np.newaxis, :]

    # generate jackknife resamples for numerator and denominator separately
    c2pt_q_jkn = generate_jackknife_resamples(c2pt_q, rwfs, bin_size)
    c2pt_0_jkn = generate_jackknife_resamples(c2pt_0, rwfs, bin_size)

    # ratio of c2pt at momenta over c2pt at momentum 0
    c2pt_ratio_jkn = c2pt_q_jkn / c2pt_0_jkn

    # nucleon mass and energy per nmax stored as gvar's
    M_gv = get_hadron_mass(ensemble, "NUCLEON")
    E_gv = compute_continuum_energy(ensemble, nmax, "NUCLEON")

    # take central values
    M = gv.mean(M_gv)
    E = gv.mean(E_gv)

    # kinematic factor
    K = 2 * E / (E + M)
    K = K[np.newaxis, np.newaxis, :, np.newaxis]
    return c2pt_ratio_jkn * K

def compute_c2pt_ratio_average(
        ensemble:           str,
        source_set:         str,
        c2pt_ratio_map_jkn: dict[int, np.ndarray],
        rwfs:               np.ndarray,
        nsquares:           list[int],
        bin_size:           int) -> dict[int, np.ndarray]:

    # average over momenta and sources
    c2pt_ratio_avg = {}
    for nmax in nsquares:
        c2pt_ratio_avg[nmax] = np.mean(
            c2pt_ratio_map_jkn[nmax], axis=(1,2)
        )
    return c2pt_ratio_avg

if __name__ == "__main__":
    ensemble   = "D251"
    source_set = "source_set1"
    nsquares   = [0, 1, 2, 3, 4, 5, 6, 8] 
    source_sink_distances = [11, 14, 16, 19]

    # threshold that cannot be exceeded, otherwise config of observable is flagged 
    SIGMA_THRESHOLD = 6

    log(f"{ensemble}: INITIALISE IDENTIFICATION OF EXCEPTIONAL CONFIGURATIONS FOR C2PT RATIOS")
    print()

    # NOTE: rwfs 
    log(f"REWEIGHTING FACTORS")
    rwfs = load_rwfs(ensemble)
    c2pt = load_c2pt(ensemble, source_set)

    print()

    log(f"IDENTIFY SMALL REWEIGHTING FACTORS")
    flag_rwf_indices = identify_small_rwf_indices(rwfs).flatten()
    N_flag_rwf       = flag_rwf_indices.shape[0]

    if N_flag_rwf > 0:
        log(f"Total {N_flag_rwf} flagged rwf(s), Config number(s) {flag_rwf_indices}")
        log(f"Exclude {N_flag_rwf} rwf(s), Config number(s) {flag_rwf_indices}")
        rwfs = np.delete(rwfs, flag_rwf_indices, axis=0)
        c2pt = np.delete(c2pt, flag_rwf_indices, axis=0)
    else:
        log("None")

    print()

    # NOTE: c2pt
    log(f"TWO-POINT FUNCTIONS")
    c2pt_ratio_map: dict[int, np.ndarray] = {}
    c2pt_ratio_avg: dict[int, np.ndarray] = {}
    log(f"Exclude config number(s) {flag_rwf_indices}")
    log("Update c2pt dimensions per nmax")

    for nmax in nsquares:
        c2pt_ratio_map[nmax] = compute_jackknife_c2pt_ratio(
            ensemble   = ensemble,
            source_set = source_set,
            c2pt       = c2pt,
            rwfs       = rwfs,
            nmax       = nmax,
            bin_size   = 1
        )

        c2pt_ratio_avg[nmax] = np.mean(c2pt_ratio_map[nmax], axis=(1,2))

    print()

    # NOTE: c2pt ratios
    log(f"TWO-POINT FUNCTION RATIOS: IDENTIFY EXCEPTIONAL CONFIGURATIONS")
    log(f"Largest source-sink distance for {ensemble}: {max(source_sink_distances)}")
    t_max = max(source_sink_distances) + 1
    log(f"Check for Euclidean timslices: {[t for t in range(t_max)]}")

    N_t = get_temporal_lattice_points(ensemble)
    exceptionals = []
    for nmax in nsquares:
        # NOTE: Check for timslices up until largest source-sink distance of c3pts 
        for t in range(19):
            exceptional = identify_exceptionals(
                c2pt_ratio_avg[nmax][:, t], SIGMA_THRESHOLD
            )
            if np.any(exceptional):
                flag_cfg_indices = np.where(exceptional)[0]
                log(f"nmax = {nmax}, t = {t}: Found {np.sum(exceptional)} exceptionals, indices = {flag_cfg_indices}")
            else: 
                continue

            exceptionals.append(exceptional)

    flagged = np.any(exceptionals, axis=0)
    N_flag_cfg = int(np.sum(flagged))
    if N_flag_cfg > 0:
        flag_cfg_indices = np.where(flagged)[0]
        log(f"Total / Indices: {N_flag_cfg} / {flag_cfg_indices}")
    else:
        log(f"No exceptionals found.")

    # exclude potential exceptionals
    if N_flag_cfg > 0:
        c2pt = np.delete(c2pt, flag_cfg_indices, axis=0)
        rwfs = np.delete(rwfs, flag_cfg_indices, axis=0)

    # Export data
    bin_sizes = [1, 2, 4, 6, 8, 10]
    for bin_size in bin_sizes:
        output_path = f"/Users/ck/projects/data/{ensemble}/c2pt_ratios/BINSIZE{bin_size}"
        os.makedirs(output_path, exist_ok=True)
        for nmax in nsquares:
            c2pt_ratio_map[nmax] = compute_jackknife_c2pt_ratio(
                ensemble   = ensemble,
                source_set = source_set,
                c2pt       = c2pt,
                rwfs       = rwfs,
                nmax       = nmax,
                bin_size   = bin_size
            )
        
            log(f"bin size = {bin_size}, nmax = {nmax}, dim = {c2pt_ratio_map[nmax].shape}")

            c2pt_ratio_avg[nmax] = np.mean(c2pt_ratio_map[nmax], axis=(1,2))
            output_file = (
                f"{ensemble}-C2PT-RATIO-JACKKNIFE"
                f"-NMAX{nmax}-BINSIZE{bin_size}"
                f"-SOURCE-SET{source_set[-1]}.pkl"
            )
            with open(f"{output_path}/{output_file}", "wb") as f:
                pickle.dump(c2pt_ratio_avg[nmax], f, protocol=pickle.HIGHEST_PROTOCOL)
                log(f"Export {output_file}.")

        print()
