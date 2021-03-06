#!/usr/bin/python env

"""
File inputs and outputs.
"""

import json
import numpy as np


def read_json(jsonfile:str) -> dict:
    """
    Read dic from json.

    Args:
        jsonfile: Json file.

    Returns:
        dict: Dictionary object.
    """
    with open(jsonfile) as f:
        dic = json.load(f)
    return dic


def get_cell_from_lammps_structure(filename:str):
    """
    Get cell from lammps structure.
    """
    from lammpkits.interfaces.lammps import get_cell_from_lammps

    lammps_input = [
            'units metal',
            'dimension 3',
            'boundary p p p',
            'atom_style atomic',
            'box tilt large',
            'read_data %s' % filename,
            ]
    cell = get_cell_from_lammps(lammps_input=lammps_input)

    return cell


def dump_cell(cell:tuple, filename:str=None, style:str='vasp'):
    """
    Dump cell into file or list for dumping yaml.

    Args:
        cell: (lattice, frac_coords, symbol).
        filename: Dump filename.
        style: Dump style. Choose from 'vasp' and 'lammps'.

    Raises:
        RuntimeError: style is not 'vasp' nor 'lammps'.
    """
    if style not in ['vasp', 'lammps']:
        raise RuntimeError("input 'style' must be 'vasp', 'lammps' or 'list'.")

    if filename is None:
        if style == 'vasp':
            fname = 'POSCAR'
        elif stype == 'lammps':
            fname = 'lammps.structure'
    else:
        fname = filename

    if style == 'vasp':
        _write_poscar(cell=cell, filename=fname)
    elif style == 'lammps':
        with open(filename, 'w') as f:
            _write_lammps_structure(
                    cell=cell,
                    filename=filename)


def _write_poscar(
        cell:tuple,
        filename:str='POSCAR'):
    """
    Write out structure to file.
    In this function, structure is not fixed
    even if its lattice basis is left handed.

    Args:
        cell: (lattice, scaled_positions, symbols).
        filename: Poscar filename.
    """
    lattice, scaled_positions, symbols = cell
    symbol_sets = list(set(symbols))
    nums = []
    idx = []
    for symbol in symbol_sets:
        index = [ i for i, s in enumerate(symbols) if s == symbol ]
        nums.append(str(len(index)))
        idx.extend(index)
    positions = np.round(np.array(scaled_positions)[idx, :],
                         decimals=9).astype(str)

    # create strings
    strings = ''
    strings += 'generated by lammpkits\n'
    strings += '1.0\n'
    for i in range(3):
        strings += ' '.join(list(np.round(
            lattice[i], decimals=9).astype(str))) + '\n'
    strings += ' '.join(symbol_sets) + '\n'
    strings += ' '.join(nums) + '\n'
    strings += 'Direct\n'
    for position in positions:
        strings += ' '.join(list(position)) + '\n'
    print("export filename:")
    print("    %s" % filename)

    with open(filename, 'w') as f:
        f.write(strings)


def _write_lammps_structure(cell:tuple,
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
