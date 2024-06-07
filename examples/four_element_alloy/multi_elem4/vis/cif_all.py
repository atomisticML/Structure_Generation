import sys
from ase.io import read,write,lammpsdata
import glob

# atomic numbers per lammps type
#Z_of_type = {1:74,2:3}
Z_of_type = {1:24,2:26,3:14,4:23}

for f in glob.glob('*.dat'):
    fprefix = f.split('.da')[0]
    atoms = lammpsdata.read_lammps_data(f, Z_of_type=Z_of_type, style="atomic")
    write('%s.cif' % fprefix,atoms)
