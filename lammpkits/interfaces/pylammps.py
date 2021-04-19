#!/usr/bin/env python

"""
This module is the interface for PyLammps.
"""

import numpy as np
from lammps import PyLammps
from pymatgen.io.lammps.data import LammpsBox


def get_cell_from_pylammps(pylmp:PyLammps, symbols:list, verbose:bool=True):
    """
    Get cell from PyLammps class object.

    Args:
        pylmp: PyLammps class object.
        symbols: Atomic symbols.
        verbose: Show detailed log.

    Returns:
        tuple: Cell.

    Raises:
        The length of symbols and atoms do not match.
    """
    assert len(symbols) == len(pylmp.atoms), \
            "The length of symbols and atoms do not match."

    box = pylmp.lmp.extract_box()
    xlo, ylo, zlo = box[0]
    xhi, yhi, zhi = box[1]
    bounds = [ [xlo, xhi], [ylo, yhi], [zlo, zhi] ]
    tilt = box[2:5]
    lattice = LammpsBox(bounds=bounds, tilt=tilt).to_lattice().matrix
    cart_coords = []
    for i in range(len(pylmp.atoms)):  # 'for atom in pylmp.atoms' not work
        cart_coords.append(pylmp.atoms[i].position)
    fixed_cart_coords = np.array(cart_coords) - np.array(box[0])
    frac_coords = np.dot(fixed_cart_coords, np.linalg.inv(lattice))

    return (lattice, frac_coords, symbols)
