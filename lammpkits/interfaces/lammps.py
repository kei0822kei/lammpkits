#!/usr/bin/env python

"""
This module is the interface for PyLammps.
"""

from typing import Union
from lammps import lammps


def get_cell_from_lammps(lammps_input:Union[list,lammps]) -> tuple:
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

    cell = get_structure_from_lammps(lammps_input=lammps_input, show_log=True)
    cell = get_cell_from_phonopy_structure(
              get_structure_from_lammps(lammps_input))

    return cell
