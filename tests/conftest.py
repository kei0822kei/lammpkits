#!/usr/bin/env python

"""
This is for pytest fixtures.
"""

import os
import pytest
from pymatgen.io.vasp.inputs import Poscar
import lammpkits
from lammpkits.interfaces.pymatgen import get_cell_from_pymatgen_structure

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(lammpkits.__file__)),
                        'tests',
                        'data')

@pytest.fixture(autouse=True, scope='session')
def mg_cell() -> tuple:
    """
    Mg hexagonal cell.

    Returns:
        tuple: Mg hexagonal cell.
    """
    pmgstructure = \
            Poscar.from_file(os.path.join(DATA_DIR, 'Mg.poscar')).structure
    cell = get_cell_from_pymatgen_structure(pmgstructure)

    return cell


@pytest.fixture(autouse=True, scope='session')
def ar_cell() -> tuple:
    """
    Ar FCC cell.

    Returns:
        tuple: Ar FCC cell.
    """
    pmgstructure = \
            Poscar.from_file(os.path.join(DATA_DIR, 'Ar.poscar')).structure
    cell = get_cell_from_pymatgen_structure(pmgstructure)

    return cell


@pytest.fixture(autouse=True, scope='session')
def mg_twinobundary_cell() -> tuple:
    """
    Mg 10-12 twinboundary cell.

    Returns:
        tuple: Mg 10-12 twinboundary cell.
    """
    pmgstructure = \
            Poscar.from_file(os.path.join(DATA_DIR, 'tb_10-12.poscar')).structure
    cell = get_cell_from_pymatgen_structure(pmgstructure)

    return cell


@pytest.fixture(autouse=True, scope='session')
def test_data_directory() -> tuple:

    return DATA_DIR
