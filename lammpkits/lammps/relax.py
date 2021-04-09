#!/usr/bin/python env

"""
Toolkits for structure relaxation with lammps.
"""

from lammps import PyLammps


class LammpsRelax():
    """
    Structure relaxation with lammps.
    """

    def __init__(
           self,
           cell:tuple,
       ):
        """
        Args:
            cell: (lattice, frac_coords, symbol)
        """
        self._lammps = PyLammps()
        self._initial_cell = cell

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

    def run_relax(self):
        """
        Run relax.
        """
