import numpy as np 

def get_momentum_list()-> list[tuple]:
    # credit: Marcel Rodekamp
    momentum_list = []

    max_nsquare: int = 8

    # loop over all positive momenta values in all three spatial directions
    for nx in range(max_nsquare+1):
        for ny in range(max_nsquare+1):
            for nz in range(max_nsquare+1):

                # ensure n^2 < n^2_max
                if nx**2 + ny**2 + nz**2 > max_nsquare: continue

                # append the positive momentum
                n = (nx,ny,nz)
                if n not in momentum_list: momentum_list.append(n)

                # go through the different combinations where a sign
                # can appear and add them to the list

                # ONE SIGN
                n = (-nx,ny,nz)
                if n not in momentum_list: momentum_list.append(n)
                n = (nx,-ny,nz)
                if n not in momentum_list: momentum_list.append(n)
                n = (nx,ny,-nz)
                if n not in momentum_list: momentum_list.append(n)

                # TWO SIGNS
                n = (nx,-ny,-nz)
                if n not in momentum_list: momentum_list.append(n)
                n = (-nx,ny,-nz)
                if n not in momentum_list: momentum_list.append(n)
                n = (-nx,-ny,nz)
                if n not in momentum_list: momentum_list.append(n)

                # THREE SIGNS
                n = (-nx,-ny,-nz)
                if n not in momentum_list: momentum_list.append(n)

    return momentum_list

def selected_momenta_indices_dict():
    mom_list = get_momentum_list()
    return {str(p): i for i, p in enumerate(mom_list)}

def get_momentum_shell(nsquare: int) -> list[tuple]:
    if nsquare == 7:
        raise ValueError(f"No data available for n2 = {nsquare}.")

    momentum_list = get_momentum_list()
    shell_dict      = selected_momenta_indices_dict()
    shell_list      = []
    shell_indices   = []

    for (nx, ny, nz) in momentum_list:
        if nx**2 + ny**2 + nz**2 == nsquare:
            shell_list.append((nx, ny, nz))
    for nvec in shell_list:
        shell_indices.append(shell_dict[str(nvec)])

    return shell_list, shell_indices

# FOR PSEUDOSCALAR AND TEMPORAL AXIAL DATA
def get_sel_momentum_list_polarisation(nsquare: int, pol: str) -> list[tuple]:
    """
    select momenta with n^2 = nsquare that can contribute to a non-zero
    pseudoscalar (or temporal axial) ratio for a given polarisation direction

    for pol = 'X','Y','Z' the corresponding component of n = (nx,ny,nz)
    has to be non-zero (since R_P^k and R_A^4,k ~ q^k)

    returns selected momenta including both signs (i.e. +nz and -nz), since the ratio
    changes sign with q^k.
    """
    if nsquare == 7:
        raise ValueError(f"No data available for nsquare = {nsquare}.")
    if nsquare < 0:
        raise ValueError("Chosen n^2 must be >= 0.")

    pol_key = pol.strip().lower()
    pol_to_idx = {"x": 0, "y": 1, "z": 2}
    if pol_key not in pol_to_idx:
        raise ValueError(f"Polarisation key must be one of 'X','Y','Z'.")

    idx = pol_to_idx[pol_key]

    # first select the momenta corresponding to nsquare = nsquare (with all sign combinations)
    shell = get_selected_momenta(nsquare)

    # then keep only those with the polarised component non-zero
    return [n for n in shell if n[idx] != 0]

def get_selected_momenta_indices(nsquare: int, pol: str) -> list[int]:
    if pol is None:
        selected_momenta = get_selected_momenta(nsquare)
    else:
        selected_momenta = get_sel_momentum_list_polarisation(nsquare, pol)

    full_momentum_list = get_momentum_list(8)
    full_momentum_indices = {n: i for i, n in enumerate(full_momentum_list)}

    selected_momentum_indices = []
    for n in selected_momenta:
        idx = full_momentum_indices.get(n)
        if idx is None:
            raise ValueError(f"Momentum {n} not found in full momentum list.")
        selected_momentum_indices.append(idx)

    return selected_momentum_indices

def get_orthogonal_momenta(nsquare: int, pol: str) -> list[int]:
    """Return list of momentum configurations (nx, ny, nz) per nsquare that are orthorgonal to polarisation direction"""
    # First get full list corresponding to  
    pol_map = {"x": 0, "y": 1, "z": 2}
    pol_idx = pol_map[pol]

    if pol not in pol_map.keys():
        raise ValueError("Polarisation must be one of x, y, z.")

    sel_mom_list     = get_selected_momenta(nsquare)
    sel_mom_dict     = selected_momenta_indices_dict(nsquare)
    orth_mom_list    = []
    orth_mom_indices = []

    for nvec in sel_mom_list:
        if nvec[pol_idx] == 0:
            orth_mom_list.append(nvec)

    for nvec in orth_mom_list:
        np.asarray(orth_mom_indices.append(sel_mom_dict[str(nvec)]))

    return orth_mom_list, orth_mom_indices
