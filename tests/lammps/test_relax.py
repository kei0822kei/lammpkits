#!/usr/bin/env python

"""
This is pytest for twinpy.structure.lattice.
"""

from lammpkits.lammps.relax import LammpsRelax
import numpy as np


def test_lammpsrelax(mg_cell):
    """
    Check Lattice.

    Todo:
        Write test for 'get_diff' and get_midpoint
        after refactering.
    """
    lmp_rlx = LammpsRelax(cell=mg_cell)
    lmp_rlx.run_relax()
    energy = lmp_rlx.energy
    np.testing.assert_allclose(energy, -3.009544466)
