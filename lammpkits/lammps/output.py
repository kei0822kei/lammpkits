#!/usr/bin/python env

"""
This module deals with lammps output.
"""

from lammpkits.interfaces.pymatgen import get_data_from_log_lammps
from lammpkits.plot.base import line_chart

class LammpsOutput():
    """
    This class deals with lammps output.
    """

    def __init__(self, logfile='log.lammps'):
        """
        Initialize.
        """
        self._logfile = logfile
        self._keys = None
        self._data = None
        self._set_data()

    def _set_data(self):
        """
        Set data of log.lammps.
        """
        keys, data = get_data_from_log_lammps(self._logfile)
        self._keys = keys
        self._data = data

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
