from ase.build import bulk,bcc100,bcc110,bcc111
from ase.io import read,write
from ase import Atoms,Atom

#atoms = bulk('MgO','rocksalt',a=4.3,cubic=True)
atoms = bulk('W','bcc',cubic=True)
#ODD numbers for symmetric bcc 100
#atoms = bcc100('Cr',size=(1,1,5))
#atoms.center(vacuum=6.0, axis=2)
big_cell = atoms * (1,1,1)

write('supercell_target.cif',big_cell)


starting_struct = atoms * (2,2,1)

write('starting.cif',starting_struct)

