#!/usr/bin/env python

"""
Interfaces for PhononLammps.
"""

import numpy as np
from phonopy import Phonopy
from phonolammps import Phonolammps
from lammpkits.file_io import dump_cell


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


def get_strings_for_phonolammps(cell:tuple,
                                filepath:str,
                                pair_style:str,
                                pair_coeff:str,
                                ):
    """
    Get lammps input for phonolammps from relax lammps input.

    Args:
        cell: Cell used for phonon calculation.
        filepath: Dump absolute file path.
        pair_style: Pair style.
        pair_coeff: Pair coefficient.
    """
    dump_cell(cell=cell, filename=filepath, style='lammps')
    symbol = cell[2][0]
    strings = [
            'units metal',
            'boundary p p p',
            'atom_style atomic',
            'box tilt large',
            'read_data %s' % filepath,
            'change_box all triclinic',
            'neigh_modify every 1 delay 0',
            'neigh_modify one 5000',
            'pair_style %s' % pair_style,
            'pair_coeff * * {} {}'.format(pair_coeff, symbol),
            ]

    return strings
