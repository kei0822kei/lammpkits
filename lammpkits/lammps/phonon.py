#!/usr/bin/python env

"""
Toolkits molecular phonon calculation with lammps.
"""

import os
import numpy as np
from lammpkits.interfaces.phonolammps import (get_phonolammps,
                                              get_phonon_from_phonolammps,
                                              get_strings_for_phonolammps)


class LammpsPhonon():
    """
    Phonon calculation with lammps.
    """

    def __init__(self,
                 cell:tuple,
                 pair_style:str,
                 pair_coeff:str,
                 dump_dir:str='.',
                 raise_dir_exists_error:bool=True,
                 ):
        """
        Init.

        Args:
            cell: Cell used for phonon calculation.
            pair_style: Pair style.
            pair_coeff: Pair coefficient.
            dump_dir: Dump directory.
            raise_dir_exists_error: Raise error if dump_dir exists.

        Raises:
            RuntimeError: Dump_dir already exists. This is activated when
                          'raise_dir_exists_error' is True.
        """
        self._cell = cell
        self._dump_dir = None
        self._set_dump_dir(dump_dir, raise_dir_exists_error)
        self._lammps_input = None
        self._set_lammps_input(cell=cell,
                               pair_style=pair_style,
                               pair_coeff=pair_coeff)
        self._phonolammps = None
        self._phonon = None

    def _set_dump_dir(self, dump_dir, raise_dir_exists_error):
        """
        Create directory and set dump dir.
        """
        if dump_dir != '.':
            if os.path.exists(dump_dir):
                if raise_dir_exists_error:
                    raise RuntimeError("directory: %s exists" % dump_dir)
            else:
                os.makedirs(dump_dir)

        self._dump_dir = dump_dir

    def _set_lammps_input(self, cell, pair_style, pair_coeff):
        """
        Set lammps input for phonon calculation.
        """
        filepath = os.path.abspath(os.path.join(self._dump_dir,
                                                'structure.lammps'))
        in_lammps = get_strings_for_phonolammps(cell=cell,
                                                filepath=filepath,
                                                pair_style=pair_style,
                                                pair_coeff=pair_coeff)

        self._lammps_input = in_lammps

    def set_phonolammps(
            self,
            supercell_matrix:np.array=np.eye(3, dtype=int),
            primitive_matrix:np.array=np.identity(3),
            show_log:bool=False,
            ):
        """
        Set phonolammps.

        Args:
            supercell_matrix: Supercell matrix.
            primitive_matrix: primitive matrix.
            show_log: If True, show log.
        """
        ph_lmp = get_phonolammps(lammps_input=self._lammps_input,
                                 supercell_matrix=supercell_matrix,
                                 primitive_matrix=primitive_matrix,
                                 show_log=show_log,
                                 )

        self._phonolammps = ph_lmp

    def run_phonon(self):
        """
        Run phonon calculation.

        Raises:
            RuntimeError: Attribute phonolammps is not set.
        """
        if self._phonolammps is None:
            raise RuntimeError("Attribute phonolammps is not set.")
        phonon = get_phonon_from_phonolammps(self._phonolammps)

        self._phonon = phonon

    def save_phonon(self, filename:str="phonopy_params.yaml"):
        """
        Save phonon.

        Args:
            filename: Save file name.

        Raises:
            RuntimeError: Attribute phonon is not set.
        """
        if self._phonon is None:
            raise RuntimeError("Attribute phonon is not set.")
        fpath = os.path.join(self._dump_dir, filename)
        self._phonon.save(filename=fpath)
