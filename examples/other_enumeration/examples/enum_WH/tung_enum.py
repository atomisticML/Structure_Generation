from ase.io import read,write
from ase import Atom,Atoms
from ase.build import bulk, fcc111, add_adsorbate
from ase.db import connect
from icet.tools import enumerate_structures
import numpy as np
import ase

def get_conc(atoms,possible_elems):
    syms = [i for i in atoms.symbols]
    concs = {elem:syms.count(elem)/len(atoms) for elem in possible_elems}
    return concs

def filter_ats(atoms):
    atlst = [at for at in atoms if at.symbol !='X']
    new_atoms = Atoms(atlst)
    new_atoms.set_cell(atoms.get_cell())
    new_atoms.set_pbc(atoms.pbc)
    return new_atoms

def change_ats(atoms,mp):
    atsyms = [s for s in atoms.symbols]
    new_syms = [mp[s] for s in atsyms]
    atoms.symbols = new_syms
    return atoms

class System_Enum:
    def __init__(self,base_species):
        self.a = None
        self.blocks = []
        self.crystal_structures = []
        self.lattice_bases = []
        self.stored_lattice = {'a':{}}
        self.base_species = base_species
        return None


    def set_crystal_structures(self,structs = ['fcc']):
        self.crystal_structures = structs

    def set_lattice_constant(self,a):
        if type(a) == float:
            self.a = [a]*len(self.crystal_structures)
        elif type(a) == list:
            self.a = a
        else:
            raise TypeError("cannot use type other than float or list for lattice constant(s)")

        return None

    def set_substitutional_blocks(self,blocks=None):
        assert len(self.crystal_structures) >= 1, "set crystal structures before defining substitutional blocks - default is 'fcc' "
        if blocks == None:
            blocks = [ [self.base_species] ] * len(self.crystal_structures)
        else:
            blocks = blocks
        self.blocks = blocks

    def enumerate_structures(self, min_int_mult, max_int_mult):
        all_structs = []
        assert self.a != None, "set_lattice_constant first, then do enumeration"
        assert self.blocks != None, "set_substitutuional_blocks first, then do structure enumeration"
        for icrystal,crystal in enumerate(self.crystal_structures):
            primitive = bulk(self.base_species[0] , crystal , a=self.a[icrystal] , cubic=False)
            print ('generating structures for crystal: %s' % crystal)
            enumerated = enumerate_structures(primitive, range(min_int_mult, max_int_mult), self.blocks[icrystal])
            sublist = [ i for i in enumerated ]
            self.stored_lattice['a'][crystal] = self.a[icrystal]
            all_structs.append(sublist)
        se.all_structs = all_structs
        return all_structs


    def compress_expand(self,axes = [0,1,2],stepsize=0.03, nsteps = 2 ):
        chem_lst = ['%s']*len(self.base_species)
        chem_str = '-'.join(b for b in chem_lst) % tuple(self.base_species)
        compressed_expanded = {icrystal: {istrct: [] for istrct in range(len(self.all_structs[icrystal])) } for icrystal in range(len(self.crystal_structures)) }
        for icrystal,crystal in enumerate(self.crystal_structures):
            unique_labelings=self.all_structs[icrystal]
            this_a = self.stored_lattice['a'][crystal]
            steps_up = np.linspace(this_a, this_a + (stepsize*nsteps), nsteps )
            steps_above = steps_up[1:]
            steps_down = np.linspace(this_a - (stepsize*nsteps), this_a, nsteps + 1 )
            steps_below = steps_down[:-1]
            all_steps = np.append(steps_below,steps_above)
            for istrct,strct in enumerate(unique_labelings):
                cell = strct.get_cell()
                scaled_pos = strct.get_scaled_positions()
                for istep in all_steps:
                    new_atoms = Atoms(strct.symbols)
                    new_cell = cell.copy()
                    for axis in axes:
                        new_cell[axis] +=  cell[axis] - istep
                    new_atoms.set_cell(new_cell)
                    new_atoms.set_scaled_positions(scaled_pos)
                    new_atoms.set_pbc(True)
                    compressed_expanded[icrystal][istrct].append(new_atoms)
                    write('ats_%s_%s_%1.3f_%d.vasp' % (chem_str,crystal,istep,istrct)  , new_atoms)

a=3.19
#base_species=[['W','Zr']],
base_species=['W','Zr'],
primitive = bulk('W' , 'bcc', a=a,cubic=False)
write('prim.cif',primitive)
enumerated = enumerate_structures(primitive, range(1,7), base_species)
pure_structs = []
for crystal in enumerated:
    atoms = change_ats(crystal,{'W':'W','Zr':'W'})
    pure_structs.append(atoms)

#comp = SymmetryEquivalenceCheck(stol=0.068)
# comp = SymmetryEquivalenceCheck(to_primitive=True)
# comp = SymmetryEquivalenceCheck()
# comp.compare(a,b)
from ase.ga.ofp_comparator import *
comp = OFPComparator()

def my_comp(a,b,check_sorted=False):
    #print (a,b, np.isclose(a.get_cell(),b.get_cell(),atol=1.e-5))
    cellflag = all([ all(a) for a in np.isclose(a.get_cell(),b.get_cell(),atol=1.e-5)])
    print ([ all(a) for a in np.isclose(a.get_cell(),b.get_cell(),atol=1.e-5)])
    try:
        if check_sorted:
            scposa = np.array(sorted(a.get_scaled_positions()))
            scposb = np.array(sorted(b.get_scaled_positions()))
            posflag = all([all(ai) for ai in np.isclose(scposa,scposb,atol=1.e-5)])
        else:
            posflag = all([all(ai) for ai in np.isclose(a.get_scaled_positions(),b.get_scaled_positions(),atol=1.e-5)])
    except:
        posflag = False
    symflag = tuple(a.symbols) == tuple(b.symbols)
    allflag = cellflag and posflag and symflag
    return allflag

def compare_wrapper(a,b):
    try:
        #res = comp.looks_like(a,b)
        res = my_comp(a,b)
    except:
        res = False
    return res

def check_against_list(atoms,lst):
    #return not any([comp.looks_like(atoms,li) for li in lst])
    #print ([compare_wrapper(atoms,li) for li in lst], (not any([compare_wrapper(atoms,li) for li in lst])))
    #return (not any([compare_wrapper(atoms,li) for li in lst]))
    return (not any([my_comp(atoms,li) for li in lst]))

keep_structs = [pure_structs[0]]
for struct in pure_structs[1:]:
    if check_against_list(struct,keep_structs):
        keep_structs.append(struct)
    else:
        pass

for icrystal,crystal in enumerate(keep_structs):
    atoms = filter_ats(crystal)
    #concs=get_conc(atoms,possible_elems=['H','W'])
    concs=get_conc(atoms,possible_elems=['W'])
    print (icrystal,concs,len(atoms))
    write('w_%03d_%02d.cif'%(icrystal,len(atoms)),atoms)
