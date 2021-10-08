#!/usr/bin/python env

"""
This module deals with lammps output.
"""

import numpy as np
from lammpkits.interfaces.pymatgen import get_data_from_log_lammps
from lammpkits.plot.base import line_chart
from lammpkits.file_io import read_json

class LammpsOutput():
    """
    This class deals with lammps output.
    """

    def __init__(self, jsonfile:str='lammpkits.json'):
        """
        Initialize.
        """
        self._lammps_input = None
        self._initial_cell = None
        self._final_cell = None
        self._keys = None
        self._data = None
        self._energy = None
        self._set_data(jsonfile)

    @property
    def initial_cell(self):
        return self._initial_cell

    @property
    def final_cell(self):
        return self._final_cell

    @property
    def keys(self):
        return self._keys

    @property
    def energy(self):
        return self._energy

    def _set_data(self, jsonfile):
        """
        Set data of log.lammps.
        """
        dump_data = read_json(jsonfile)
        self._lammps_input = dump_data['lammps_input']
        self._keys = dump_data['log_lammps']['keys']
        self._data = np.array(dump_data['log_lammps']['data'])
        self._initial_cell = dump_data['initial_cell']
        self._final_cell = dump_data['initial_cell']
        for cell in [self._initial_cell, self._final_cell]:
            cell[0] = np.array(cell[0])
            cell[1] = np.array(cell[1])
        self._energy = self.get_data_from_key('PotEng')[-1]

    def get_keys(self):
        """
        Get keys.
        """
        return self._keys

    def get_data(self):
        """
        Get data.
        """
        return self._data

    def get_data_from_key(self, key):
        """
        Get data from key.
        """
        if not key in self._keys:
            raise RuntimeError("key: {} is not exist.".format(key))
        idx = self._keys.index(key)
        key_data = self._data[idx]

        return key_data

    def plot_transition(self, ax, key:str):
        """
        Plot transition per steps. xlabel is always steps.

        Args:
            ax: Ax of matplotlib.
            key: Key for ylabel.
        """
        xlabel = 'Step'
        ylabel = key
        xdata = self.get_data_from_key(xlabel)
        ydata = self.get_data_from_key(ylabel)
        line_chart(ax=ax,
                   xdata=xdata,
                   ydata=ydata,
                   xlabel=xlabel,
                   ylabel=ylabel,
                   )
        return ax
