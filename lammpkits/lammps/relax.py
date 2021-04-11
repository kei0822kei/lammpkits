#!/usr/bin/python env

"""
Toolkits for structure relaxation with lammps.
"""

import os
import tempfile
import numpy as np
from lammps import PyLammps
import lammpkits
from lammpkits.file_io import write_lammps_structure
from lammpkits.interfaces.pylammps import get_cell_from_pylammps


class LammpsRelax():
    """
    Structure relaxation with lammps.
    """

    def __init__(self):
        self._lammps = PyLammps()
        self._initial_cell = None
        self._lammps_input = []
        self._lammps_potential_symbols = None
        self._is_run_finished = False

    def _check_run_is_finished(self):
        """
        Check lammps run is finished.

        Raises:
            RuntimeError: Lammps run is not finished.
        """
        if not self._is_run_finished:
            raise RuntimeError("Lammps run is not finished.")

    def _check_run_is_not_finished(self):
        """
        Check lammps run not is finished.

        Raises:
            RuntimeError: Lammps run is finished.
        """
        if self._is_run_finished:
            raise RuntimeError("Lammps run is finished.")

    @property
    def lammps(self):
        """
        PyLammps class object.
        """
        return self._lammps

    def get_initial_cell(self):
        """
        Get initial cell.

        Returns:
            tuple: Initial cell.
        """
        return self._initial_cell

    @property
    def lammps_input(self):
        """
        Lammps input.
        """
        return self._lammps_input

    def print_lammps_input(self):
        """
        Print lammps input.
        """
        print('\n'.join(self._lammps_input))

    def add_structure(self, cell:tuple):
        """
        Add structure to lammps input.

        Args:
            cell: (lattice, frac_coords, symbol)
        """
        self._check_run_is_not_finished()
        structure_file = _dump_cell(cell)
        strings = [
                'units metal',
                'dimension 3',
                'boundary p p p',
                'atom_style atomic',
                'atom_modify map array',
                'box tilt large',
                'read_data %s' % structure_file,
                ]
        self._initial_cell = cell
        self._lammps_potential_symbols = tuple(set(cell[2]))
        self._lammps_input.extend(strings)

    def add_potential_from_string(self, pair_style:str, pair_coeff:str):
        """
        Add potential to lammps input.

        Args:
            pair_style: Key 'pair_style' setting for lammps input.
            pair_coeff: Key 'pair_coeff' setting for lammps input.
        """
        self._check_run_is_not_finished()
        strings = [ pair_style, pair_coeff ]
        self._lammps_input.extend(strings)

    def add_potential_from_database(self, pair_style:str, pot_file:str):
        """
        Add potential to lammps input from database.

        Args:
            pair_style: Key 'pair_style' setting for lammps input.
            pot_file: Potential file path from potentials directory.

        Notes:
            Potential file is read from
            /path/to/lammpkits/potentials/dir/pot_file.
            You can create potentials directory by
            ```
            mkdir /path/to/lammpkits/potentials
            ```
            and add potentials.
        """
        self._check_run_is_not_finished()
        pot_dir = os.path.join(os.path.dirname(os.path.dirname(
                                               lammpkits.__file__)),
                               'potentials')
        potential_file_path = os.path.join(pot_dir, pot_file)
        pair_coeff = 'pair_coeff * * {} {}'.format(
                potential_file_path,
                ' '.join(self._lammps_potential_symbols),
                )
        self.add_potential_from_string(
                pair_style=pair_style,
                pair_coeff=pair_coeff,
                )

    def add_relax_settings(self):
        """
        Add relax settings.
        """
        self._check_run_is_not_finished()
        strings = [
                'neighbor 0.3 bin',
                'fix 1 all box/relax aniso 0.0 couple xy vmax 0.01',
                'variable energy equal pe',
                'minimize 0.0 1.0e-8 1000 100000',
                ]
        self._lammps_input.extend(strings)

    def run_lammps(self, use_mpi:bool=False):
        """
        Run lammps.

        Args:
            use_mpi: Whether use mpi or not.

        Notes:
            When use_mpi=True, the script must be run with mpirun like
            ```
            mpirun -np 4 python script.py
            ```
        """
        self._check_run_is_not_finished()
        fname = _dump_strings(self._lammps_input)
        if use_mpi:
            from mpi4py import MPI
        self._lammps.file(fname)
        if use_mpi:
            MPI.Finalize()
        self._is_run_finished = True

    def get_final_cell(self):
        """
        Get final cell

        Returns:
            tuple: Final cell.
        """
        self._check_run_is_finished()
        return get_cell_from_pylammps(
                pylmp=self._lammps,
                symbols=self._initial_cell[2],
                )

    def get_energy(self):
        """
        Get energy after relax.

        Args:
            float: Final energy.
        """
        self._check_run_is_finished()
        return self._lammps.variables['energy'].value

    def get_forces(self):
        """
        Get forces acting on atoms after relax.

        Args:
            np.array: Forces acting on atoms.
        """
        atoms = self._lammps.atoms
        forces = []
        for i in range(len(atoms)):  # 'for atom in atoms' does not work
            forces.append(atoms[i].force)
        return np.array(forces)

    def get_lammps_input_for_phonolammps(self):
        """
        Get lammps input for phonolammps.

        Args:
            list: List of lammps commands for phonoLAMMPS.
        """
        self._check_run_is_finished()
        structure_file = _dump_cell(self.get_final_cell())
        strings = [
                'units metal',
                'boundary p p p',
                'atom_style atomic',
                'box tilt large',
                'read_data %s' % structure_file,
                ]
        pot_strings = [ s for s in self._lammps_input
                            if 'pair_style' in s or 'pair_coeff' in s ]
        strings.extend(pot_strings)
        strings.append('neighbor 0.3 bin')

        return strings


def _dump_cell(_cell):
    """
    Return template file path.
    """
    _, tmp_fname = tempfile.mkstemp()
    write_lammps_structure(
            cell=_cell,
            filename=tmp_fname)
    return tmp_fname


def _dump_strings(strings):
    """
    Return template file path.
    """
    _, tmp_fname = tempfile.mkstemp()
    with open(tmp_fname, 'w') as f:
        f.write('\n'.join(strings))
    return tmp_fname
