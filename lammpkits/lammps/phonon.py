#!/usr/bin/python env

"""
Toolkits molecular phonon calculation with lammps.
"""


class LammpsPhonon():
    """
    Phonon calculation with lammps.
    """

    def __init__(self,
                 in_lammps:str,
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
        self._lammps_input = in_lammps

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

