#!/usr/bin/env python

import numpy as np
from pymatgen.io.vasp.inputs import Poscar
from phonopy import Phonopy
from phonolammps import Phonolammps
from lammpkits.lammps.static import LammpsStatic
from lammpkits.interfaces.pymatgen import get_cell_from_pymatgen_structure


def main():
    pmgstructure = Poscar.from_file("tb_10-12_orig.poscar").structure
    cell = get_cell_from_pymatgen_structure(pmgstructure)
    mg_pair_style = 'pair_style mlip_pair'
    pot_file = "Mg/mlp/pot_icsd/pair-1/mlp.lammps"

    lmp_rlx = LammpsStatic()
    lmp_rlx.add_structure(cell=cell)
    lmp_rlx.add_potential_from_database(
            pair_style=mg_pair_style,
            pot_file=pot_file,
            )
    lmp_rlx.add_variables()
    lmp_rlx.add_relax_settings(is_relax_lattice=True)
    lmp_rlx.run_lammps()

    ####  phonon
    strings = lmp_rlx.get_lammps_input_for_phonolammps()
    supercell_matrix = np.array([[3, 0, 0], [0, 1, 0], [0, 0, 1]])
    primitive_matrix = np.eye(3)

    # from lammpkits.lammps.static import _dump_strings
    # fpath = _dump_strings(strings)
    # print(fpath)
    # # print("hoge"*1000)

    phl = Phonolammps(strings,
                      supercell_matrix=supercell_matrix,
                      primitive_matrix=primitive_matrix,
                      show_log=False,
                      show_progress=False)
    unitcell = phl.get_unitcell()
    force_constants = phl.get_force_constants()
    phonon = phl.get_phonopy_phonon()
    phl.plot_phonon_dispersion_bands()

    # It seems necessary to add some codes to phonolammps to get phonon class object
    # which is already set forces.
    # Can generate FORCE_SETS from force_constants easily?
    phonon = Phonopy(unitcell=unitcell, supercell_matrix=supercell_matrix,
                     primitive_matrix=primitive_matrix,
                     )  # Probably it is not enough if you need FORCE_SETS

if __name__ == '__main__':
    main()
