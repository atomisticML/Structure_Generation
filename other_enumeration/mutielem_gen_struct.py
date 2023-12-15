import ase
from ase.io import read,write
from ase import Atoms,Atom
from ase.neighborlist import primitive_neighbor_list
from ase.data import atomic_numbers
from ase.ga.utilities import closest_distances_generator
from ase.data import atomic_numbers
from ase.build import bulk
from ase.spacegroup import crystal
from hnf import *
import math

def get_num_dense(atoms):
	cell_i = atoms.get_cell()
	vol = np.dot(np.cross(cell_i[0],cell_i[1]) , cell_i[2])
	ndense = len(atoms)/vol
	return ndense,vol

#s=Atoms('', pbc=True)

chems = ['Ta','V']
desired_comps = {'Ta':0.9, 'V':0.1}
#HNF matrix trace (unique supercell from hermite-normal-form)
N = 60
full_hnf_mats = get_hnfs(hnf_trs = [N])
hnf_mats_sub = limit_mats_len(full_hnf_mats,N,tol =0.16 )
hnf_mats = hnf_mats_sub.copy()
#np.random.seed(23331)
cell_id = np.random.choice(range(len(hnf_mats)))
cell_mult = hnf_mats[cell_id]
print ('this hnf',cell_mult)
a=3.4

tmpats = bulk('Ta','bcc',a=a,cubic=False)
write('crystal.cif',tmpats)

dns,prim_vol = get_num_dense(tmpats)


cell_c = tmpats.get_cell().copy()
cell_col = cell_c.T
hnf_cell = np.matmul(cell_col,cell_mult).T

target_vol = np.dot(np.cross(hnf_cell[0],hnf_cell[1]),hnf_cell[2])
mult = target_vol/prim_vol
desired_size = int(mult)

print( hnf_cell,desired_size)


def flip_one_atom(atoms,types):
	new_atoms = atoms.copy()
	flip_ind = np.random.randint(0,len(atoms))
	flip_current = new_atoms[flip_ind].symbol
	excluded = [typ for typ in types if typ != flip_current]
	flip_to_ind = np.random.randint(0,len(excluded))
	flip_to_type = excluded[flip_to_ind]
	new_atoms[flip_ind].symbol = flip_to_type
	return new_atoms

def flip_N_atoms(atoms,types,fraction=0.25):
	pert_inds = np.random.choice(range(len(atoms)),size=int(len(atoms)*fraction) )
	new_atoms = atoms.copy()
	for pert_ind in pert_inds:
		flip_ind = np.random.randint(0,len(atoms))
		flip_current = new_atoms[flip_ind].symbol
		excluded = [typ for typ in types if typ != flip_current]
		flip_to_ind = np.random.randint(0,len(excluded))
		flip_to_type = excluded[flip_to_ind]
		new_atoms[flip_ind].symbol = flip_to_type
	return new_atoms

def add_atom(atoms,symbols,tol = 0.5):
	blmin = closest_distances_generator(atom_numbers=[atomic_numbers[symbol] for symbol in symbols] + [atomic_numbers['Ne']], ratio_of_covalent_radii=0.5)
	def readd():
		symbol = np.random.choice(symbols)
		rnd_pos_scale = np.random.rand(1,3)
		rnd_pos = np.matmul(atoms.get_cell(),rnd_pos_scale.T)
		rnd_pos = rnd_pos.T[0]
		new_atom = Atom('Ne',rnd_pos)
		tst_atoms = atoms.copy()
		tst_atoms.append(new_atom)
		tst_atoms.wrap()
		rc = 5.
		
		atinds = [atom.index for atom in tst_atoms]
		at_dists = {i:[] for i in atinds}
		all_dists = []
		nl = primitive_neighbor_list('ijdD',pbc=tst_atoms.pbc,positions=tst_atoms.positions ,cell=atoms.get_cell(),cutoff=rc)
		bond_types = {i:[] for i in atinds}
		for i,j in zip(nl[0],nl[-1]):
			at_dists[i].append(j)
		for i,j in zip(nl[0],nl[1]):
			bond_types[i].append( (atomic_numbers[tst_atoms[i].symbol] , atomic_numbers[tst_atoms[j].symbol])  )
		return symbol, tst_atoms, at_dists, rnd_pos, bond_types
	symbol, tst_atoms , at_dists , rnd_pos, bond_types = readd()
	bondtyplst = list(bond_types.keys())
	syms = [tst_atom.symbol for tst_atom in tst_atoms]
	tst_id = syms.index('Ne')
	tst_dists = at_dists[tst_id]
	tst_bonds = bond_types[tst_id]
	conds = all([ np.linalg.norm(tst_dist) >=  blmin[(atomic_numbers[symbol] , tst_bonds[i][1])] for i,tst_dist in enumerate(tst_dists)])
	while not conds:
		symbol , tst_atoms, at_dists , rnd_pos, bond_types = readd()
		syms = [tst_atom.symbol for tst_atom in tst_atoms]
		tst_id = syms.index('Ne')
		tst_dists = at_dists[tst_id]
		tst_bonds = bond_types[tst_id]
		#conds = all([ np.linalg.norm(tst_dist) >= tol for tst_dist in tst_dists])
		#conds = all([ np.linalg.norm(tst_dist) >= blmin[tst_bonds[i]] for i,tst_dist in enumerate(tst_dists)])
		print (tst_dists)
		print (tst_bonds)
		conds = all([ np.linalg.norm(tst_dist) >=  blmin[(atomic_numbers[symbol] , tst_bonds[i][1])] for i,tst_dist in enumerate(tst_dists)])
	atoms.append(Atom(symbol,rnd_pos))
	return atoms


def get_comp(atoms,symbols):
	comps = {symbol: 0.0 for symbol in symbols}
	counts = {symbol: 0 for symbol in symbols}
	atsymbols = [atom.symbol for atom in atoms]
	for atsymbol in atsymbols:
		counts[atsymbol] +=1
	for symbol in symbols:
		comps[symbol] = counts[symbol]/len(atoms)
	return comps

initial_a = Atoms('',pbc=True)
initial_a.set_cell(hnf_cell)

current_comps = {'Ta':0.0,'V':0.0}
current_size = 0

import time

comp_conds = all([current_comps[chem] == desired_comps[chem] for chem in chems])
these_atoms = initial_a.copy()
max_iter = 1000
this_iter = 0 
while current_size <= desired_size or not comp_conds and this_iter <= max_iter:
	t = 1000 * time.time() # current time in milliseconds
	np.random.seed(int(t) % 2**32)
	if current_size < desired_size:
		tst_ats = add_atom(these_atoms,chems)
		tst_comps =  get_comp(tst_ats, chems)
		comp_conds = all([tst_comps[chem] == desired_comps[chem] for chem in chems])
		these_atoms = tst_ats.copy()
		current_size = len(these_atoms)
		Qi = [np.abs(tst_comps[chem] - desired_comps[chem]) for chem in chems]
		Qi = round(np.sum(Qi),8)
	elif current_size == desired_size:
		#tst_ats = flip_one_atom(these_atoms,chems)
		tst_ats = flip_N_atoms(these_atoms,chems,fraction=np.random.rand())
		tst_comps =  get_comp(tst_ats, chems)
		comp_conds = all([tst_comps[chem] == desired_comps[chem] for chem in chems])
		tstQi = [np.abs(tst_comps[chem] - desired_comps[chem]) for chem in chems]
		tstQi = round(np.sum(tstQi),8)
		#print (len(these_atoms),tst_comps,tstQi)
		if tstQi <= Qi:
			these_atoms = tst_ats.copy()
			Qi = tstQi
			print ('in composition loop',tst_comps,comp_conds)
	if current_size==desired_size and comp_conds:
		break
	this_iter += 1
print (these_atoms)

write('initialatoms.cif',these_atoms)
