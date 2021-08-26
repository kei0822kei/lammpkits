#!/usr/bin/env python

"""
Interfaces for PhononLammps.
"""

import numpy as np
from phonolammps import Phonolammps


def get_phonolammps(lammps_input:list,
                    supercell_matrix:np.array=np.eye(3, dtype=int),
                    primitive_matrix:np.array=np.identity(3),
                    show_log:bool=True,
                    ):
    """
    Docstring for get_phonolammps.

    Args:
        supercell_matrix: Supercell matrix.
        primitive_matrix: primitive matrix.
        show_log: If True, show log.
    """
    ph_lmp = Phonolammps(lammps_input=lammps_input,
                         supercell_matrix=supercell_matrix,
                         primitive_matrix=primitive_matrix,
                         show_log=show_log,
                         symmetrize=True,  # symmetrize is True by default
                         )

    return ph_lmp


def get_phonon_from_phonolammps(phonolammps) -> Phonopy:
    """
    Get Phonopy class object from PhonoLammps.

    Args:
        phonolammps: Phonlammps class object.
    """
    unitcell = phonolammps.get_unitcell()
    force_constants = phonolammps.get_force_constants()
    supercell_matrix = phonolammps.get_supercell_matrix()
    primitive_matrix = phonolammps.get_primitve_matrix()
    phonon = Phonopy(unitcell=unitcell,
                     primitive_matrix=primitive_matrix,
                     supercell_matrix=supercell_matrix)
    dataset = phonolammps.get_force_constants(include_data_set=True)[1]
    phonon.dataset = dataset
    phonon.produce_force_constants()

    return phonon
