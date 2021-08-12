#!/usr/bin/python env

"""
Toolkits molecular phonon calculation with lammps.
"""


class LammpsPhonon():
    """
    Phonon calculation with lammps.
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
