#!/usr/bin/env python

"""
This is pytest for twinpy.structure.lattice.
"""

from pymatgen.io.vasp.inputs import Poscar
from lammpkits.interfaces.pymatgen import get_cell_from_pymatgen_structure
from lammpkits.lammps.static import LammpsStatic

ar_pair_style = 'pair_style lj/cut 10'
ar_pair_coeff = 'pair_coeff * * 0.01042 3.4 10'

pmgstructure = \
        Poscar.from_file('Ar.poscar').structure
ar_cell = get_cell_from_pymatgen_structure(pmgstructure)

lmp_rlx = LammpsStatic()
lmp_rlx.add_structure(cell=ar_cell)
lmp_rlx.add_potential_from_string(
        pair_style=ar_pair_style,
        pair_coeff=ar_pair_coeff,
        )
lmp_rlx.add_thermo(thermo=100)
lmp_rlx.add_variables(add_energy=True, add_stress=True)
lmp_rlx.add_relax_settings(is_relax_lattice=True)
lmp_rlx.run_lammps()
energy = lmp_rlx.get_energy()
print("energy: %f" % energy)
