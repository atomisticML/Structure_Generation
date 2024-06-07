import sys
from ase.io import read,write,lammpsdata

atoms = lammpsdata.read_lammps_data(sys.argv[1], Z_of_type={1:24}, style="atomic")
#atoms = read(sys.argv[1],format='lammps-data')
write('tmp.cif',atoms)
