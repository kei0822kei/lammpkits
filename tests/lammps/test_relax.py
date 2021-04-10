#!/usr/bin/env python

"""
This is pytest for twinpy.structure.lattice.
"""

from lammpkits.lammps.relax import LammpsRelax
import numpy as np


def test_lammpsrelax(ar_cell):
    """
    Check Lattice.
    """
    pair_style = 'pair_style lj/cut 10'
    pair_coeff = 'pair_coeff * * 0.01042 3.4 10'

    lmp_rlx = LammpsRelax()
    lmp_rlx.add_structure(cell=ar_cell)
    lmp_rlx.add_potential_from_string(
            pair_style=pair_style,
            pair_coeff=pair_coeff,
            )
    lmp_rlx.add_relax_settings()
    lmp_rlx.run_lammps()
    energy = lmp_rlx.energy
    np.testing.assert_allclose(energy, -0.346153, atol=1e-5)
