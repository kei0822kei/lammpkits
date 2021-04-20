#!/usr/bin/python env

"""
Toolkits molecular static calculation for lammps.
"""

import os
import tempfile
from pprint import pprint
import numpy as np
from lammps import lammps
import lammpkits
from lammpkits.file_io import write_lammps_structure
from lammpkits.interfaces.lammps import get_cell_from_lammps


class LammpsStatic():
    """
    Molecular static calculation for lammps.
    """

    def __init__(self):
        self._lammps = lammps()
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

    @property
    def is_run_finished(self):
        """
        Is lammps run finished.
        """
        return self._is_run_finished

    def get_initial_cell(self) -> tuple:
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
            cell: (lattice, frac_coords, symbol).
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
                'change_box all triclinic',
                'neigh_modify every 1 delay 0',
                'neigh_modify one 5000',
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

    def add_thermo(self, thermo:int=100):
        """
        Add thermo settings.

        Args:
            thermo: Int of thermo setting.
        """
        self._check_run_is_not_finished()
        strings.append('thermo %d' % thermo)
        strings.append('thermo_style custom step pe press '
                       + 'pxx pyy pzz pxy pxz pyz lx ly lz vol')
        strings.append('thermo_modify norm no')  # Do not normalize

    def add_variables(self,
                      add_energy:bool=True,
                      add_stress:bool=True,
                      ):
        """
        Add variables.

        Args:
            add_energy: If True, variable energy is set.
            add_stress: If True, variable stress is set.

        Todo:
            Currently add_stress is not working.
        """
        self._check_run_is_not_finished()
        strings = []

        if add_energy:
            strings.append('variable energy equal pe')

        if add_stress:
            strings.append('variable pxx0 equal pxx')
            strings.append('variable pyy0 equal pyy')
            strings.append('variable pzz0 equal pzz')
            strings.append('variable pyz0 equal pyz')
            strings.append('variable pxz0 equal pxz')
            strings.append('variable pxy0 equal pxy')

        self._lammps_input.extend(strings)

    def add_relax_settings(self, is_relax_lattice:bool=True):
        """
        Add relax settings.

        Args:
            is_relax_lattice: If True, lattice is relaxed.
        """
        self._check_run_is_not_finished()
        strings = []
        if is_relax_lattice:
            strings.append('fix 1 all box/relax tri 0')
        strings.append('min_style cg')
        strings.append('minimize 0.0 1.0e-8 100000000 1000000000')
        self._lammps_input.extend(strings)

    def run_lammps(self, verbose:bool=True):
        """
        Run lammps.

        Args:
            verbose: If True, show detailed information.
        """
        self._check_run_is_not_finished()
        fname = _dump_strings(self._lammps_input)
        if verbose:
            print("Dump lammps input to %s" % fname)
            print("Run lammps with the inputs:")
            pprint(self._lammps_input)
        self._lammps.file(fname)
        self._is_run_finished = True

    def get_final_cell(self) -> tuple:
        """
        Get final cell

        Returns:
            tuple: Final cell.
        """
        self._check_run_is_finished()
        cell = get_cell_from_lammps(self._lammps)

        return cell

    def get_energy(self):
        """
        Get energy after relax.

        Args:
            float: Final energy.
        """
        self._check_run_is_finished()
        return self._lammps.extract_variable('energy', None, 0)

    def get_stress(self) -> list:
        """
        Get stress after relax.

        Args:
            list: List of [ Pxx Pyy Pzz Pyz Pzx Pxy ]

        Todo:
            Future return 3x3 numpy array.
        """
        self._check_run_is_finished()
        keys = ['pxx0', 'pyy0', 'pzz0', 'pyz0', 'pxz0', 'pxy0']
        stress = [ self._lammps.extract_variable(key, None, 0)
                       for key in keys ]
        return stress

    def get_forces(self) -> np.array:
        """
        Get forces acting on atoms after relax.

        Args:
            np.array: Forces acting on atoms.
        """
        self._check_run_is_finished()
        n = self._lammps.get_natoms()
        f = self._lammps.extract_atom('f', 3)
        forces = np.array([[f[i][0],f[i][1],f[i][2]] for i in range(n)])

        return forces

    def get_lammps_input_for_phonolammps(self) -> list:
        """
        Get lammps input for phonolammps.

        Args:
            list: List of lammps commands for phonoLAMMPS.
        """
        self._check_run_is_finished()
        structure_file = _dump_cell(self.get_final_cell())
        pot_strings = [ s for s in self._lammps_input
                            if 'pair_style' in s or 'pair_coeff' in s ]

        strings = [
                'units metal',
                'boundary p p p',
                'atom_style atomic',
                'box tilt large',
                'read_data %s' % structure_file,
                'change_box all triclinic',
                'neigh_modify every 1 delay 0',
                'neigh_modify one 5000',
                ]
        strings.extend(pot_strings)

        return strings


def _dump_cell(cell:tuple) -> str:
    """
    Return template file path.

    Args:
        cell: (lattice, frac_coords, symbol).

    Returns:
        str: Dumped file path.
    """
    _, tmp_fname = tempfile.mkstemp()
    write_lammps_structure(
            cell=cell,
            filename=tmp_fname)

    return tmp_fname


def _dump_strings(strings:list) -> str:
    """
    Return template file path.

    Args:
        strings: List of strings.

    Returns:
        str: Dumped file path.
    """
    _, tmp_fname = tempfile.mkstemp()
    with open(tmp_fname, 'w') as f:
        f.write('\n'.join(strings))

    return tmp_fname
