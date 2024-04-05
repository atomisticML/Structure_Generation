import sys
from ase.io import read,write,lammpsdata
import glob

for f in glob.glob('*.dat'):
    fprefix = f.split('.da')[0]
    atoms = lammpsdata.read_lammps_data(f, Z_of_type={1:24}, style="atomic")
    write('%s.cif' % fprefix,atoms)
