#!/usr/bin/python env

"""
File inputs and outputs.
"""

def write_lammps_structure(cell:tuple,
                           filename:str):
    """
    Write lammps structure.

    Args:
        cell: (lattice, frac_coords, symbols).
        filename: Output filename.
    """
    from lammpkits.interfaces.pymatgen import get_pymatgen_structure
    from pymatgen.io.lammps.data import LammpsData

    pmgstruct = get_pymatgen_structure(cell=cell)
    lmp_data = LammpsData.from_structure(pmgstruct, atom_style='atomic')
    lmp_data.write_file(filename=filename)
