#!/usr/bin/env python

import os
import argparse
import numpy as np
from pymatgen.io.vasp.inputs import Poscar
from lammpkits.api_lammpkits import Lammpkits
from lammpkits.interfaces.pymatgen import get_cell_from_pymatgen_structure

def get_argparse():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--posfile',
                        type=str,
                        default='POSCAR',
                        help="POSCAR file.")
    parser.add_argument('--pair_style',
                        type=str,
                        help="choose 'mlip_pair' or 'mlip_coeff'.")
    parser.add_argument('--pair_coeff',
                        type=str,
                        help="mlp.lammps file path")
    parser.add_argument('--supercell_matrix',
                        type=str,
                        help="Supercell maxtrix. ex. '3 1 1'.")
    argments = parser.parse_args()
    return argments


def main(
        posfile,
        pair_style,
        pair_coeff,
        supercell_matrix,
        ):

    # cell
    pmgstructure = Poscar.from_file(posfile).structure
    cell = get_cell_from_pymatgen_structure(pmgstructure)

    # lammps
    lmpkits = Lammpkits(cell=cell)
    lmpkits.run_lammps_phonon(
            cell=cell,
            pair_style=pair_style,
            pair_coeff=pair_coeff,
            supercell_matrix=supercell_matrix,
            )


if __name__ == '__main__':
    args = get_argparse()
    supercell_matrix = np.diag(list(map(int, args.supercell_matrix.split())))
    main(posfile=args.posfile,
         pair_style=args.pair_style,
         pair_coeff=args.pair_coeff,
         supercell_matrix=supercell_matrix,
         )
