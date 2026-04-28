from collections import namedtuple

Lattice = namedtuple('Lattice', ['n_space', 'n_time'])

_LATTICE_DIMENSIONS = {
    # A
    # B
    # C
    'C101': Lattice(48,  96),
    'C102': Lattice(48,  96),
    # D
    'D251': Lattice(64, 128),
    # E
    # H
    'H106': Lattice(32,  96),
    'H107': Lattice(32,  96),
    # J
    # N
    'N201': Lattice(48, 128),
    'N204': Lattice(48, 128),
    'N304': Lattice(48, 128),
    # S
    # U
    # X
}