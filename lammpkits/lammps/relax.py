#!/usr/bin/python env

"""
Toolkits for structure relaxation with lammps.
"""

import tempfile
from lammps import PyLammps
from lammpkits.file_io import write_lammps_structure


class LammpsRelax():
    """
    Structure relaxation with lammps.
    """

    def __init__(self):
        self._lammps = PyLammps()
        self._initial_cell = None
        self._final_cell = None
        self._energy = None
        self._lammps_input = []

    @property
    def lammps(self):
        """
        PyLammps class object.
        """
        return self._lammps

    @property
    def initial_cell(self):
        """
        Initial cell.
        """
        return self._initial_cell

    @property
    def energy(self):
        """
        Initial cell.
        """
        return self._energy

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
        for line in self._lammps_input:
            print(line)

    def add_structure(self, cell:tuple):
        """
        Add structure to lammps input.

        Args:
            cell: (lattice, frac_coords, symbol)
        """
        def _dump_cell(_cell):
            """
            Return template file path.
            """
            _, tmp_fname = tempfile.mkstemp()
            write_lammps_structure(
                    cell=_cell,
                    filename=tmp_fname)
            return tmp_fname

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
        self._lammps_input.extend(strings)

    def add_potential_from_string(self, pair_style:str, pair_coeff:str):
        """
        Add potential to lammps input.

        Args:
            pair_style: Key 'pair_style' setting for lammps input.
            pair_coeff: Key 'pair_coeff' setting for lammps input.
        """
        strings = [ pair_style, pair_coeff ]
        self._lammps_input.extend(strings)

    def add_relax_settings(self):
        """
        Add relax settings.
        """
        strings = [
                'neighbor 0.3 bin',
                'fix 1 all box/relax aniso 0.0 couple xy vmax 0.01',
                'variable energy equal pe',
                'minimize 0.0 1.0e-8 1000 100000',
                ]
        self._lammps_input.extend(strings)

    def run_lammps(self):
        """
        Run lammps.

        Raises:
            RuntimeError: Structure relaxation is alread conducted.
        """
        def _set_vals(lmp):
            """
            store results from lammps.
            """
            self._energy = lmp.variables['energy'].value

        if self._energy:
            raise RuntimeError("Structure relaxation is already conducted.")

        for line in self._lammps_input:
            self._lammps.command(line)

        _set_vals(self._lammps)
