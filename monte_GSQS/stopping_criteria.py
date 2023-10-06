from __future__ import print_function
import sys,os
import ctypes
import math
import numpy as np
from ase.io import read,write
from ase import Atoms,Atom
from ase.ga.utilities import closest_distances_generator
from ase.data import atomic_numbers
from ase.neighborlist import primitive_neighbor_list
from lammps import lammps, LMP_TYPE_ARRAY, LMP_STYLE_GLOBAL
from lammps import PyLammps
# get MPI settings from LAMMPS
lmp = lammps()
L = PyLammps(ptr=lmp)


start = 'fcc_4.cif'
file_prefix = 'iter_%d' % 0
atoms = read(start)*(2,2,1)

write('%s.data' % file_prefix,atoms,format='lammps-data')

fname = file_prefix
L.command("clear")
L.command("info all out log")
L.command('units  metal')
L.command('atom_style  atomic')
L.command("boundary   p p p")
L.command("atom_modify	map hash")
L.command('neighbor  2.3 bin')
# boundary
L.command('boundary  p p p')
# read atoms
L.command('read_data  %s.data' % fname )
L.command('mass  1 180.94788')

# potential settings
L.command("pair_style	zero 5.7")
L.command("pair_coeff	* *")

# define compute pace

L.command("compute   pace all pace coupling_coefficients.yace 1 0")

# run
#L.command("fix	  1 all nve")
#L.command("timestep     0.00025")
L.command("thermo	    1")
L.run(0)
#L.run(1, "pre no post no")
#L.command("run	    0")

L_pace = lmp.numpy.extract_compute("pace", LMP_STYLE_GLOBAL, LMP_TYPE_ARRAY)
descriptors_start_decomp = L_pace[ : len(atoms), : -1]
#descriptors_start = np.sum(descriptors_start_decomp,axis=0)
descriptors_start = np.average(descriptors_start_decomp,axis=0)
descriptors_last = descriptors_start

natoms = len(atoms)
def check_atoms(tst_id,tst_atoms,symbols=['Ta']):
	blmin = closest_distances_generator(atom_numbers=[atomic_numbers[symbol] for symbol in symbols], ratio_of_covalent_radii=0.35)
	def readd():
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
		return at_dists, bond_types
	at_dists , bond_types = readd()
	bondtyplst = list(bond_types.keys())
	syms = [tst_atom.symbol for tst_atom in tst_atoms]
	tst_dists = at_dists[tst_id]
	tst_bonds = bond_types[tst_id]
	conds = all([ np.linalg.norm(tst_dist) >=  blmin[(atomic_numbers[tst_atoms[tst_id].symbol] , tst_bonds[i][1])] for i,tst_dist in enumerate(tst_dists)])
	return conds

def get_atoms(positions,symbols,cell):
	atoms_inter = Atoms()
	for atiind in range(len(positions)):
		pos = positions[atiind]
		atoms_inter.append(Atom(symbols[atiind],pos))
	atoms_inter.set_cell(cell)
	atoms_inter.set_pbc(True)
	return atoms_inter

def run_lmp(niterations,deltamove,beta,stopping_criteria=0.001,last_residual=0.0):
	naccept = 0
	residuals = []
	descriptors_target = np.load('target_descriptors.npy')
	descriptors_target = descriptors_target#*len(read('supercell.cif'))
	last_abs_diffs = [ np.abs(ii-tt)  for ii,tt in zip(descriptors_start,descriptors_target)]
	last_residual = np.sum(last_abs_diffs)
	last_residual=last_residual
	for iiter in range(niterations):
		iatom = np.random.choice(range(len(atoms)))
		current_atom = L.atoms[iatom]

		x0, y0, z0 = current_atom.position

		dx = deltamove * np.random.uniform(-1, 1)
		dy = deltamove * np.random.uniform(-1, 1)
		dz = deltamove * np.random.uniform(-1, 1)

		current_atom.position = (x0+dx, y0+dy, z0+dz)
		#print ([x0,y0,z0],current_atom.position)
		L.run(1, "pre no post no")
		#L.command("run	    0")

		L_pace = lmp.numpy.extract_compute("pace", LMP_STYLE_GLOBAL, LMP_TYPE_ARRAY)
		descriptors_tst_decomp = L_pace[ : len(atoms), : -1]
		#descriptors_tst = np.sum(descriptors_tst_decomp,axis=0)
		descriptors_tst = np.average(descriptors_tst_decomp,axis=0)

		abs_diffs = [np.abs(ii - kk) for ii,kk in zip(descriptors_tst, descriptors_target)]
		abs_diffs = np.array(abs_diffs)
		
		#descriptors.append(descriptors_tst)

		#tst_residual = np.sum(abs_diffs) 
		tst_residual = np.sum(abs_diffs) 
		resid_diff =  (tst_residual - last_residual)

		full_tst_residual =  (abs_diffs - last_abs_diffs)

		boltzmann_dists = np.exp(-(full_tst_residual)*beta)

		boltzmann_dists = list(boltzmann_dists)
		many_rands = list(np.random.rand(len(boltzmann_dists)))

		meet_ratio = 0.75
		discard_ratio = 1 - meet_ratio

		rndi = np.random.rand()

		tst_positions = []
		for atiind in range(len(atoms)):
			ati = L.atoms[atiind]
			tst_positions.append(ati.position)
		atoms_tst = get_atoms(tst_positions,['Ta']*len(tst_positions),atoms.get_cell())
		logical3 = check_atoms(tst_id=iatom,tst_atoms=atoms_tst)

		logical1 = tst_residual <= last_residual
		bolt1 = np.exp(-resid_diff*beta)
		if np.isinf(bolt1):
			bolt1 = 0.0
		logical2 = rndi <= bolt1
		#logical3 = True #np.abs(np.log10(tst_residual/last_residual)) <= 6
		#boltzmann_flags = [ r <= bd for r,bd in zip(many_rands,boltzmann_dists) ]
		#some_flag = boltzmann_flags.count(True) >= int(len(boltzmann_flags)/1.25)
		#some_flag2 = boltzmann_flags.count(True) >= int(len(boltzmann_flags)/(1+discard_ratio))
		#simple_lt = [abs_diffsi <= last_absdiffsi for abs_diffsi, last_absdiffsi in zip(abs_diffs,last_abs_diffs)]
		#some_flag1 = simple_lt.count(True) >= int(len(boltzmann_flags)/(1+(4*discard_ratio)))
		if logical1 and logical3:#tst_residual <= last_residual:
		#if some_flag1 and tst_residual <= last_residual:
			naccept += 1
			#print (iiter,tst_residual,last_residual)
			last_residual = tst_residual
			last_abs_diffs = abs_diffs
		#elif some_flag2:
		elif logical2 and logical3 and not logical1:# np.random.rand() <= np.exp((resid_diff)*beta):
			#print (resid_diff,rndi,bolt1,logical2)
			#print (logical1)
			naccept += 1
			last_residual = tst_residual
			last_abs_diffs = abs_diffs
		else:
			#print ('not accepted')
			current_atom.position = (x0, y0, z0)
		print (beta,iiter,last_residual, naccept, bolt1,np.std(residuals[-200:]))
		if iiter >=10:
			if np.std(residuals[-200:]) <= stopping_criteria:
				break
		residuals.append(last_residual)
	these_positions = []
	for atiind in range(len(atoms)):
		ati = L.atoms[atiind]
		these_positions.append(ati.position)
	atoms_inter = get_atoms(these_positions,['Ta']*len(these_positions),atoms.get_cell())
	write('beta_%f_atoms.cif' % beta,atoms_inter)
	return residuals,last_abs_diffs


#thermodynamic betas
#betas = [5.e-2,5.e-1,5.e0,5.e1,5.e2,5.e3]
betas = [5.e-2,5.e-1,5.e0,5.e1,5.e2,5.e2]
#N iterations 
#iters = [100,200,500,500,1000]
#iters = [100]*6
iters = [300]*2 + [500]*2 + [3000] +  [2000]
#delta moves
deltamoves = [0.5, 0.5, 0.125, 0.075, 0.05, 0.005]

for ibeta,beta in enumerate(betas):
	try:
		last_residual = np.load('residuals%s.npy' % str(ibeta))[0]
	except FileNotFoundError:
		last_residual = np.sum(descriptors_start)
	resids, last_abs_diffs = run_lmp(iters[ibeta],deltamoves[ibeta],beta,last_residual=last_residual)
	np.save('residuals%s.npy' % str(ibeta+1),np.array(resids))
	print (last_abs_diffs)

