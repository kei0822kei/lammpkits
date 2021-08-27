#!/usr/bin/python env

"""
Toolkits molecular static calculation for lammps.
"""

import os
from pprint import pprint
import json
import numpy as np
from lammps import lammps
import lammpkits
from lammpkits.file_io import dump_cell
from lammpkits.interfaces.lammps import get_cell_from_lammps
from lammpkits.interfaces.pymatgen import get_data_from_log_lammps


class LammpsStatic():
    """
    Molecular static calculation for lammps.
    """

    def __init__(self,
                 dump_dir:str='.',
                 raise_dir_exists_error:bool=True,
                 ):
        """
        Init.

        Args:
            dump_dir: Dump directory.
            raise_dir_exists_error: Raise error if dump_dir exists.

        Raises:
            RuntimeError: Dump_dir already exists. This is activated when
                          'raise_dir_exists_error' is True.
        """
        self._dump_dir = None
        self._set_dump_dir(dump_dir, raise_dir_exists_error)
        self._lammps = lammps(
                cmdargs=[
                    '-log',
                    os.path.join(dump_dir, 'log.lammps'),
                    ]
                )
        self._initial_cell = None
        self._lammps_input = []
        self._lammps_potential_symbols = None
        self._is_run_finished = False
        self._lattice_type = None

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

    def add_string(self, string):
        self._lammps_input.append(string)

    def add_structure(self,
                      cell:tuple,
                      dump_filename:str='initial_structure.lammps'):
        """
        Add structure to lammps input.

        Args:
            cell: (lattice, frac_coords, symbol).
        """
        self._check_run_is_not_finished()
        structure_fpath = os.path.join(
                os.getcwd(),
                self._dump_dir,
                dump_filename)
        dump_cell(cell=cell, filename=structure_fpath, style='lammps')
        strings = [
                'units metal',
                'dimension 3',
                'boundary p p p',
                'atom_style atomic',
                'atom_modify map array',
                'box tilt large',
                'read_data %s' % structure_fpath,
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
        symbol = self._initial_cell[2][0]
        pot_dir = os.path.join(os.path.dirname(os.path.dirname(
                                               lammpkits.__file__)),
                               'potentials')
        potential_file_path = os.path.join(pot_dir, symbol, pot_file)
        pair_coeff = 'pair_coeff * * {} {}'.format(
                potential_file_path,
                ' '.join(self._lammps_potential_symbols),
                )
        self.add_potential_from_string(
                pair_style=pair_style,
                pair_coeff=pair_coeff,
                )

    def add_thermo(self, thermo:int=1000):
        """
        Add thermo settings.

        Args:
            thermo: Int of thermo setting.
        """
        self._check_run_is_not_finished()
        strings = []

        strings.append('thermo %d' % thermo)
        strings.append('thermo_style custom step pe press '
                       + 'pxx pyy pzz pxy pxz pyz lx ly lz vol fmax fnorm')
        strings.append('thermo_modify norm no')  # Do not normalize

        self._lammps_input.extend(strings)

    def add_dump(self, dump_steps:int=10, basedir:str='cfg'):
        """
        Dump cells.

        Args:
            dump_steps: Dump cells every input value steps.
            basedir: Base directory for storing cells.
        """
        self._check_run_is_not_finished()
        if self._initial_cell is None:
            raise RuntimeError("Structure is not set.")
        dump_structure_dir = os.path.join(self._dump_dir, basedir)

        strings = []
        strings.append(
                "shell mkdir %s" % dump_structure_dir)
        dump = "d1 all cfg {} {}/run.*.cfg mass type xs ys zs id type".format(
                   dump_steps, dump_structure_dir)
        strings.append('dump %s' % dump)
        self._lammps_input.extend(strings)

    def add_fix_box_relax(self, keyvals:dict):
        """
        Add fix box/relax command. Available keys and values are shown in
        https://docs.lammps.org/fix_box_relax.html.

        Args:
            keyvals: If keyvals = {'x': 0, 'y': 0}, then fix box/relax command
                     is added as 'fix box/relax x 0 y 0'.
        """
        keys_vals = []
        for key in keyvals:
            keys_vals.extend([key, keyvals[key]])
        string = 'fix f1 all box/relax %s' % ' '.join(map(str, keys_vals))
        self._lammps_input.append(string)

    def add_unfix(self):
        """
        Add unfix.
        """
        self._check_run_is_not_finished()
        self._lammps_input.extend('unfix f1')

    def add_minimize(self,
                     etol:float=1e-10,
                     ftol:float=1e-10,
                     maxiter:int=100000,
                     maxeval:int=100000,
                     ):
        """
        Add relax settings.

        Note:
            The details of input args are shown in
            https://docs.lammps.org/minimize.html.
        """
        self._check_run_is_not_finished()
        strings = []
        strings.append('reset_timestep 0')
        strings.append('min_style cg')
        strings.append('minimize {} {} {} {}'.format(
            etol, ftol, maxiter, maxeval))
        self._lammps_input.extend(strings)

    def run_lammps(self, verbose:bool=True, lammps_filename:str='in.lammps'):
        """
        Run lammps.

        Args:
            verbose: If True, show detailed information.
        """
        self._check_run_is_not_finished()
        lammps_fpath = os.path.join(
                os.getcwd(),
                self._dump_dir,
                lammps_filename)
        with open(lammps_fpath, 'w') as f:
            f.write('\n'.join(self._lammps_input))
        if verbose:
            print("Dump lammps input to %s" % lammps_fpath)
            print("Run lammps with the inputs:")
            pprint(self._lammps_input)
        self._lammps.file(lammps_fpath)
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
        return self._lammps.get_thermo('pe')

    def get_stress(self) -> list:
        """
        Get stress after relax.

        Args:
            list: List of [ Pxx Pyy Pzz Pyz Pzx Pxy ]

        Todo:
            Future return 3x3 numpy array.
        """
        self._check_run_is_finished()
        keys = ['pxx', 'pyy', 'pzz', 'pyz', 'pxz', 'pxy']
        stress = [ self._lammps.get_thermo(key) for key in keys ]
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

    # def get_lammps_input_for_phonolammps(
    #         self,
    #         dump_filename:str='final_structure.lammps') -> list:
    #     """
    #     Get lammps input for phonolammps.

    #     Args:
    #         list: List of lammps commands for phonoLAMMPS.
    #     """
    #     strings = get_lammps_input_for_phonolammps()
    #     self._check_run_is_finished()
    #     structure_fpath = os.path.join(os.getcwd(),
    #                                    self._dump_dir,
    #                                    dump_filename)
    #     final_cell = self.get_final_cell()
    #     dump_cell(cell=final_cell, filename=structure_fpath, style='lammps')
    #     pot_strings = [ s for s in self._lammps_input
    #                         if 'pair_style' in s or 'pair_coeff' in s ]

    #     strings = [
    #             'units metal',
    #             'boundary p p p',
    #             'atom_style atomic',
    #             'box tilt large',
    #             'read_data %s' % structure_fpath,
    #             'change_box all triclinic',
    #             'neigh_modify every 1 delay 0',
    #             'neigh_modify one 5000',
    #             ]
    #     strings.extend(pot_strings)

    #     return strings

    def as_dict(self):
        """
        Get dict including all lammps inputs and outputs.
        """
        self._check_run_is_finished()
        logfile = os.path.join(self._dump_dir, 'log.lammps')
        keys, data = get_data_from_log_lammps(logfile)
        log_lammps = {'keys': keys, 'data': data}

        dic = {
            'lammps_input': self._lammps_input,
            'initial_cell': self._initial_cell,
            'final_cell': self.get_final_cell(),
            'log_lammps': log_lammps,
        }

        return dic

    def dump_lammps(self, filename='lammpkits.json'):
        """
        Dump lammps.

        Args:
            filename: Dump file name.
        """

        dic = self.as_dict()
        dic['initial_cell'] = [ dic['initial_cell'][0].tolist(),
                                dic['initial_cell'][1].tolist(),
                                dic['initial_cell'][2] ]
        dic['final_cell'] = [ dic['final_cell'][0].tolist(),
                                dic['final_cell'][1].tolist(),
                                dic['final_cell'][2] ]
        dic['log_lammps']['data'] = dic['log_lammps']['data'].tolist()
        with open(filename, 'w') as f:
            json.dump(dic, f, indent=4)
