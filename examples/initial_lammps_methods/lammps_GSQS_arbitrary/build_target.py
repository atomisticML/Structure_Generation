from ase.build import bulk,bcc100,bcc110,bcc111
from ase.io import read,write
from ase import Atoms,Atom

atoms = bulk('W','bcc',cubic=True)
atoms[1].symbol = 'Be'
write('supercell_target.cif',atoms)
