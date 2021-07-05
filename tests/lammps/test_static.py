#!/usr/bin/env python

"""
This is pytest for twinpy.structure.lattice.
"""

from lammpkits.lammps.static import LammpsStatic
import numpy as np


def test_lammpsrelax(ar_cell, mg_cell):
    """
    Check Lattice.
    """
    # Ar
    ar_pair_style = 'pair_style lj/cut 10'
    ar_pair_coeff = 'pair_coeff * * 0.01042 3.4 10'

    lmp_rlx = LammpsStatic(dump_dir='./dump',
                           raise_dir_exists_error=False)
    lmp_rlx.add_structure(
            cell=ar_cell,
            dump_filename='initial_structure.lammps',
            )
    lmp_rlx.add_potential_from_string(
            pair_style=ar_pair_style,
            pair_coeff=ar_pair_coeff,
            )
    lmp_rlx.add_thermo()
    lmp_rlx.add_variables()
    lmp_rlx.add_relax_settings(is_relax_lattice=True)
    lmp_rlx.print_lammps_input()
    lmp_rlx.run_lammps(lammps_filename='in.lammps')
    energy = lmp_rlx.get_energy()
    np.testing.assert_allclose(energy, -0.346153, atol=1e-5)

    # # Mg
    # from lammpkits.file_io import write_poscar
    # mg_pair_style = 'pair_style mlip_pair'
    # pot_file = "Mg/mlp/pot/pair-8/mlp.lammps"

    # lmp_rlx = LammpsStatic()
    # lmp_rlx.add_structure(cell=mg_cell)
    # lmp_rlx.add_potential_from_database(
    #         pair_style=mg_pair_style,
    #         pot_file=pot_file,
    #         )
    # lmp_rlx.add_relax_settings()
    # lmp_rlx.run_lammps()
    # energy = lmp_rlx.get_energy()
    # np.testing.assert_allclose(energy, -2.985058, atol=1e-5)
    # write_poscar(lmp_rlx.get_initial_cell(), "/home/mizo/initial.poscar")
    # write_poscar(lmp_rlx.get_final_cell(), "/home/mizo/final.poscar")
