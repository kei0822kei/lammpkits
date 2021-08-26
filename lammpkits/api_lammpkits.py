#!/usr/bin/python env

"""
This module provides API for lammpkits.
"""

import os
from lammpkits.file_io import dump_cell
from lammpkits.lammps.static import LammpsStatic
from lammpkits.lammps.output import LammpsOutput

class Lammpkits():
    """
    API for lammpkits.
    """

    def __init__(self, cell:tuple):
        """
        Initialize.
        """
        self._initial_cell = cell
        self._final_cell = None
        self._statics = None
        self._outputs = None

    @property
    def initial_cell(self):
        """
        Initial cell.
        """
        return self._initial_cell

    @property
    def final_cell(self):
        """
        Final cell.
        """
        return self._final_cell

    @property
    def statics(self):
        """
        List of LammpsStatic class objects.
        """
        return self._statics

    @property
    def outputs(self):
        """
        List of LammpsOutput class objects.
        """
        return self._outputs

    def run_static_calc(
            self,
            pair_style:str,
            pot_file:str,
            minimize_settings:list,
            base_dir:str='.',
            dump_steps:int=50,
            ):
        """
        Set LammpsStatic and LammpsOutput.

        Args:
            pair_style: Key 'pair_style' setting for lammps input.
            pot_file: Potential file path from potentials directory.
            minimize_settings: Minimize settings. See Examples.
            base_dir: Dump base directory.
            dump_steps: Dump every 'dump_steps'.
            is_dump_lammps: If True, dump lammps results.

        Examples:
            Here is the example for 'minimize_settings'.

            >>> minimize_settings = [
                            {
                                'box_relax': {'aniso': 0, 'couple': 'xy'},
                                'minimize': {'etol': 1e-6,
                                             'ftol': 1e-6,
                                             'maxiter': 10000,
                                             'maxeval': 10000}
                            },
                            {
                                'box_relax': {'tri': 0},
                                'minimize': {'etol': 1e-10,
                                             'ftol': 1e-10,
                                             'maxiter': 10000,
                                             'maxeval': 10000}
                            },
                            {
                                'minimize': {'etol': 1e-10,
                                             'ftol': 1e-10,
                                             'maxiter': 10000,
                                             'maxeval': 10000}
                            }
                        ]
        """
        statics = []
        outputs = []
        initial_cell = self._initial_cell
        for i, lmp_args in enumerate(minimize_settings):
            dump_dir = os.path.join(base_dir, 'minimize_'+str(i))
            lmp_stc = LammpsStatic(dump_dir=dump_dir)
            lmp_stc.add_structure(cell=initial_cell)
            dump_cell(cell=lmp_stc.get_initial_cell(),
                      filename=os.path.join(dump_dir, 'initial_cell.poscar'))
            lmp_stc.add_potential_from_database(pair_style=pair_style, pot_file=pot_file)
            lmp_stc.add_thermo(thermo=dump_steps)
            lmp_stc.add_dump(dump_steps=dump_steps, basedir='cfg')
            if 'box_relax' in lmp_args.keys():
                lmp_stc.add_fix_box_relax(keyvals=lmp_args['box_relax'])
            lmp_stc.add_minimize(**lmp_args['minimize'])
            lmp_stc.run_lammps()
            lmp_stc.dump_lammps(
                    filename=os.path.join(dump_dir,'lammpkits.json'))
            final_cell = lmp_stc.get_final_cell()
            dump_cell(cell=final_cell,
                      filename=os.path.join(dump_dir, 'final_cell.poscar'))
            logfile = os.path.join(dump_dir, 'log.lammps')
            lmp_out = LammpsOutput(logfile)
            statics.append(lmp_stc)
            outputs.append(lmp_out)
            initial_cell = final_cell

        self._final_cell = final_cell
        self._statics = statics
        self._outputs = outputs

    def run_lammps_phonon(
            self,
            ):
        """
        Run lammps phonon.
        """
