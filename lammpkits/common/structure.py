#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This module deals with hexagonal twin structure.
"""

from pprint import pprint
import numpy as np
from phonopy.structure.atoms import atom_data, symbol_map
from twinpy.properties.hexagonal import (get_hcp_atom_positions,
                                         check_cell_is_hcp)
from twinpy.properties.twinmode import TwinIndices
from twinpy.structure.lattice import CrystalLattice


def get_numbers_from_symbols(symbols:list):
    """
    Get atomic numbers from symbols.

    Args:
        symbols: Atomic symbols.
    """
    numbers = [ symbol_map[symbol] for symbol in symbols ]
    return numbers


def get_symbols_from_numbers(numbers:list):
    """
    Get symbols from atomic numbers.

    Args:
        numbers: Atomic numbers.
    """
    symbols = [ atom_data[number][1] for number in numbers ]
    return symbols


def check_same_cells(first_cell:tuple,
                     second_cell:tuple,
                     raise_error:bool=False,
                     atol:float=1e-6) -> bool:
    """
    Check first cell and second cell are same.

    Args:
        first_cell: First cell.
        second_cell: Second cell.
        raise_error: If True, raise assrtion error.

    Returns:
        bool: Return True if two cells are same.
    """
    is_lattice_same = np.allclose(first_cell[0], second_cell[0], atol=atol)
    is_scaled_positions_same = np.allclose(
            first_cell[1], second_cell[1], atol=atol)
    is_symbols_same = (first_cell[2] == second_cell[2])
    is_same = (is_lattice_same
               and is_scaled_positions_same
               and is_symbols_same)
    if not is_same and raise_error:
        np.testing.assert_allclose(first_cell[0], second_cell[0], atol=atol)
        np.testing.assert_allclose(first_cell[1], second_cell[1], atol=atol)
        assert (first_cell[2] == second_cell[2])

    return is_same


def get_atom_positions_from_lattice_points(lattice_points:np.array,
                                           atoms_from_lp:np.array) -> np.array:
    """
    Get atom positions by embedding primitive atoms to lattice points.
    Both lattice points and atom positions must be cartesian coordinates.

    Args:
        lattice_points: Lattice points.
        atoms_from_lp: Atoms from lattice_points.

    Returns:
        np.array: atom positions
    """
    scaled_positions = []
    for lattice_point in lattice_points:
        atoms = lattice_point + atoms_from_lp
        scaled_positions.extend(atoms.tolist())
    return np.array(scaled_positions)
