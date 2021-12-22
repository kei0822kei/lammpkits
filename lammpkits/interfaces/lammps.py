#!/usr/bin/env python

"""
This module is the interface for PyLammps.
"""

from typing import Union
import numpy as np
from lammps import lammps


def get_cell_from_lammps(lammps_input:Union[list,lammps],
                         move_atoms_into_unitcell:bool=True) -> tuple:
    """
    Get cell from lammps input object or strings.

    Args:
        lammps_input: This is allowed two styles. One is list which has
                      lammps_input strings necessary for extracting structure.
                      The other is lammps class object after structure is set.

    Returns:
        tuple: Extracted cell.
    """
    from phonolammps.iofile import get_structure_from_lammps
    from lammpkits.interfaces.phonopy import get_cell_from_phonopy_structure

    lattice, scaled_positions, symbols = \
      get_cell_from_phonopy_structure(get_structure_from_lammps(lammps_input))
    scaled_positions = np.round(scaled_positions, decimals=8)
    if move_atoms_into_unitcell:
        scaled_positions %= 1.

    return (lattice, scaled_positions, symbols)
