from ase.build import bulk,bcc100,bcc110,bcc111
from ase.io import read,write
from ase import Atoms,Atom
import numpy as np

def get_target_comp_s(desired_size,desired_comps):
    num_ats = {key:int(desired_comps[key]*desired_size) for key in list(desired_comps.keys())}
    symbols =[]
    for key,nrepeat in num_ats.items():
        symbols.extend([key]*nrepeat)
    np.random.shuffle(symbols)
    return symbols

desired_comps = {'Cr':0.7,'Fe':0.2, 'Si':0.05,'V':0.05}

atoms = bulk('Cr','bcc',cubic=True)
#atoms = atoms*(12,12,12)
m = 24 #48*2
atoms = atoms*(m,m,m)
desired_size = len(atoms)
new_comps = {elem:int(round(desired_size*cm))/desired_size for elem,cm in desired_comps.items()}
#new_comps = {elem:int(desired_size*cm)/desired_size for elem,cm in desired_comps.items()}
print(desired_comps)
print(new_comps)
syms = get_target_comp_s(desired_size,new_comps)
print('desired size vs syms',desired_size,len(syms))
atoms.symbols = syms

#write('supercell_target.cif',atoms)
write('supercell_target.xyz',atoms)
