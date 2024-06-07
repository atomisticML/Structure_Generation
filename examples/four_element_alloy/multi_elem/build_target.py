from ase.build import bulk,bcc100,bcc110,bcc111
from ase.io import read,write
from ase import Atoms,Atom

atoms = bulk('Cr','bcc',cubic=True)
atoms = atoms*(2,1,1)
atoms[1].symbol = 'Fe'
atoms[2].symbol = 'Si'
atoms[3].symbol = 'V'

write('supercell_target.cif',atoms)
