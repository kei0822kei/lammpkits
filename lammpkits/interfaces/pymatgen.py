#!/usr/bin/env python

"""
This module is the interface for pymatgen.
"""

from pymatgen.core.structure import Structure


def get_pymatgen_structure(cell:tuple) -> Structure:
    """
    Get pymatgen structure from cell.

    Args:
        cell: Cell (lattice, scaled_positions, symbols).
    """
    return Structure(lattice=cell[0],
                     coords=cell[1],
                     species=cell[2])


def get_cell_from_pymatgen_structure(pmgstructure:Structure) -> tuple:
    """
    Get cell from pymatgen.

    Args:
        pmgstructure: Pymatgen structure.

    Returns:
        tuple: Cell (lattice, scaled_positions, symbols).
    """
    lattice = pmgstructure.lattice.matrix
    scaled_positions = pmgstructure.frac_coords
    symbols = [ specie.value for specie in pmgstructure.species ]

    return (lattice, scaled_positions, symbols)


def get_data_from_log_lammps(logfile='log.lammps') -> tuple:
    """
    Get data from log.lammps.

    Args:
        logfile: Lammps logfile.

    Return:
        tuple: (keys [list], data [np.array]).
    """
    from pymatgen.io.lammps.outputs import parse_lammps_log
    df = parse_lammps_log(logfile)[0]
    keys = df.keys().tolist()
    data = df.T.to_numpy()

    return (keys, data)
