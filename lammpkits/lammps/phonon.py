#!/usr/bin/python env

"""
Toolkits molecular phonon calculation with lammps.
"""

from lammpkits.interfaces.phonolammps import (get_phonolammps,
                                              get_phonon_from_phonolammps)


class LammpsPhonon():
    """
    Phonon calculation with lammps.
    """

    def __init__(self,
                 in_lammps:list,
                 dump_dir:str='.',
                 raise_dir_exists_error:bool=True,
                 ):
        """
        Init.

        Args:
            in_lammps: Lammps input for phonon.
            dump_dir: Dump directory.
            raise_dir_exists_error: Raise error if dump_dir exists.

        Raises:
            RuntimeError: Dump_dir already exists. This is activated when
                          'raise_dir_exists_error' is True.
        """
        self._dump_dir = None
        self._set_dump_dir(dump_dir, raise_dir_exists_error)
        self._lammps_input = in_lammps
        self._phonolammps = None
        self._phonon = None

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

    def set_phonolammps(
            self,
            supercell_matrix:np.array=np.eye(3, dtype=int),
            primitive_matrix:np.array=np.identity(3),
            show_log:bool=True,
            ):
        """
        Set phonolammps.

        Args:
            supercell_matrix: Supercell matrix.
            primitive_matrix: primitive matrix.
            show_log: If True, show log.
        """
        ph_lmp = get_phonolammps(lammps_input=self._lammps_input,
                                 supercell_matrix=supercell_matrix,
                                 primitive_matrix=primitive_matrix)

        self._phonolammps = ph_lmp

    def run_phonon(self):
        """
        Run phonon calculation.

        Raises:
            RuntimeError: Attribute phonolammps is not set.
        """
        if self._phonolammps is None:
            raise RuntimeError("Attribute phonolammps is not set.")
        phonon = get_phonon_from_phonolammps(self._phonolammps)

        self._phonon = phonon

    def save_phonon(self, filename:str="phonopy_params.yaml"):
        """
        Save phonon.

        Args:
            filename: Save file name.

        Raises:
            RuntimeError: Attribute phonon is not set.
        """
        if self._phonon is None:
            raise RuntimeError("Attribute phonon is not set.")
        self._phonon.save(filename=filename)
