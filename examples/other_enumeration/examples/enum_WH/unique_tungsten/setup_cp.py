from ase.io import read,write,iread
from qe_tools import *
import os

""" &system
    ibrav=0, nat=2, ntyp=1,nspin=2,
    ecutwfc  =  60.000000000000,
    ecutrho  =  480,
    starting_magnetization(1)  =  5.000000000000D-03,
    nr1b=12, nr2b=12, nr3b=24,
 /
 &electrons
    electron_dynamics='cg', electron_damping=0.6,
    startingwfc='random', ampre=0.1,
    emass=700., emass_cutoff=3.,
    electron_temperature='nose',
"""
def setup_cp(atoms,prefix='test'):
    qei = QEInput(atoms)
    qei.default_settings()
    fname = '%s.in' % prefix
    cpstr1 = """ &control
    calculation='cp',
    restart_mode='from_scratch',
    verbosity='high',
    nstep=%d, iprint=%d, isave=%d,
    dt=1.0,
    saverho=.false.,
    disk_io='low',
    outdir='./out',"""
    cpstr2="""&system
    ibrav=0, nat=%d, ntyp=%d,nspin=2,
    ecutwfc  =  60.000000000000,
    ecutrho  =  720,
    nr1b=24, nr2b=24, nr3b=24,"""
    cpstr3="""
&electrons
    electron_dynamics='cg', electron_damping=0.6,
    startingwfc='random', ampre=0.1,
    emass=700., emass_cutoff=3.,
    electron_temperature='nose',
/
&ions
    ion_dynamics='damp', ion_temperature='nose',
    ion_damping=0.15,
    tempw=300,
/
&cell
  cell_dynamics='damp-pr',
  cell_damping=10.15,
  cell_dofree='all',
/

    """
    qei.write_cp_input(cpstr1,cpstr2,cpstr3,fname)

import sys
#atomsitr = [at for at in iread(sys.argv[1])]
#atoms_final = atomsitr[-1]
#setup_cp(atoms_final)
#fname = sys.argv[1]
for ifn,fname in enumerate(sorted(glob.glob('w_*.cif'))):
    atoms = read(fname)
    prefix = fname.split('.')[0] + '_md'
    dirname = 'run_%03d' % ifn
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
    os.chdir(dirname)
    setup_cp(atoms,prefix)
    os.chdir('../')
        
