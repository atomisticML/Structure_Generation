from crystal_enum import *
from ase.io import read,write

atoms = read('supercell_target.xyz')

def get_symmetries(atoms):
    this_lattice_typ = get_cell_type(atoms,parent_only=False)
    parent_lattice_typ = get_cell_type(atoms,parent_only=True)
    return this_lattice_typ, parent_lattice_typ

syms = get_symmetries(atoms)
print(syms)
