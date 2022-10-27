
from __future__ import print_function
import sys,os
import ctypes
import numpy as np
from ase.io import read,write
from subprocess import call
from lammps import lammps, LMP_TYPE_ARRAY, LMP_STYLE_GLOBAL

cleanup=False

# get MPI settings from LAMMPS
def run_struct(in_atoms,fname,type_map={1:'Ca',2:'Mg'},masses={1: 40.078,2:24.305}):
	map_type = {value:key for key,value in type_map.items()}
	lmp = lammps()
	me = lmp.extract_setting("world_rank")
	nprocs = lmp.extract_setting("world_size")
	species = sorted(list(set([a.symbol for a in in_atoms])))
	species_types = [map_type[spec] for spec in species]
	inp_masses = [masses[typ] for typ in species_types]
	cmds = ["-screen", "none", "-log", "none"]
	lmp = lammps(cmdargs = cmds)

	#def set_atoms(atoms,atid=0):
	#	write('iter_%d.data' % atid,atoms,format='lammps-data')
	#	lmp.command('read_data  iter_%d.data' % atid )
	#	lmp.command('mass  1 180.94788')
	#	lmp.command(f"run {nsteps}")

	def run_lammps(dgradflag):

		# simulation settings
		fname = file_prefix
		lmp.command("clear")
		lmp.command("info all out log")
		lmp.command('units  metal')
		lmp.command('atom_style  atomic')
		#lmp.command("atom_modify	map hash")
		lmp.command('neighbor  2.3 bin')
		# boundary
		#lmp.command('boundary  p p p')
		# read atoms
		lmp.command('read_data  %s.data' % fname )
		for ityp,species_typ in enumerate(species_types):
			imass = inp_masses[ityp]
			lmp.command('mass  %d %f' % (ityp+1,imass))

		# potential settings

		lmp.command(f"pair_style 	zero 4.5")
		lmp.command(f"pair_coeff 	* *")

		# define compute pace

		if dgradflag:
			lmp.command(f"compute 	desc all pace coupling_coefficients.yace 1 1 ")
		else:
			lmp.command(f"compute 	desc all pace coupling_coefficients.yace 1 0")

		# run

		lmp.command(f"thermo 		1")
		lmp.command(f"thermo_style    custom step temp pe ke etotal c_desc[1][1]")
		#lmp.command(f"run {nsteps}")
		lmp.command(f"run  0")

	nsteps = 0
	# declare compute pace variables
	bikflag = 1
	dgradflag = 0

	run_lammps(dgradflag)

	lmp_pace = lmp.numpy.extract_compute("desc", LMP_STYLE_GLOBAL, LMP_TYPE_ARRAY)
	L_pace = lmp_pace.copy()
	lmp.close()
	del lmp
	return L_pace

def rescale_num(num,a=0.5,b=1.0):
	minval=0.0
	maxval=1.0
	prefac = b-a
	ratio = (num - minval)/(maxval-minval)
	rescaled = (prefac*ratio) + a
	return rescaled

import sys

def cost_func_1(l_pace,descriptors_target,atoms):
	descriptors_tst_decomp = l_pace[ : len(atoms), : -1]
	#descriptors_tst = np.sum(descriptors_tst_decomp,axis=0)
	descriptors_tst = np.average(descriptors_tst_decomp,axis=0)
	abs_diffs = [np.abs(ii - kk) for ii,kk in zip(descriptors_tst, descriptors_target)]
	abs_diffs = np.array(abs_diffs)
	tst_residual = np.sum(abs_diffs)
	return tst_residual,abs_diffs

def cost_func_2(l_pace,descriptors_target,atoms):
	descriptors_tst_decomp = l_pace[ : len(atoms), : -1]
	descriptors_tst = np.var(descriptors_tst_decomp,axis=0)
	abs_diffs = [np.abs(ii - kk) for ii,kk in zip(descriptors_tst, descriptors_target)]
	abs_diffs = np.array(abs_diffs)
	tst_residual = np.sum(abs_diffs)
	return tst_residual,abs_diffs

#f = sys.argv[-1]

def get_comp(atoms,symbols):
	comps = {symbol: 0.0 for symbol in symbols}
	counts = {symbol: 0 for symbol in symbols}
	atsymbols = [atom.symbol for atom in atoms]
	for atsymbol in atsymbols:
		counts[atsymbol] +=1
	for symbol in symbols:
		comps[symbol] = counts[symbol]/len(atoms)
	return comps

def perturb_one_atom(atoms,scale=0.5):
	new_atoms = atoms.copy()
	pert_ind = np.random.randint(0,len(atoms))
	perturbation = np.random.rand(1,3)[0]
	posneg = 2.*(perturbation - np.min(perturbation))/np.ptp(perturbation)-1
	posneg *= scale
	new_atoms[pert_ind].x += posneg[0]
	new_atoms[pert_ind].y += posneg[1]
	new_atoms[pert_ind].z += posneg[2]
	return new_atoms


def flip_one_atom(atoms,types):
	new_atoms = atoms.copy()
	flip_ind = np.random.randint(0,len(atoms))
	flip_current = new_atoms[flip_ind].symbol
	excluded = [typ for typ in types if typ != flip_current]
	flip_to_ind = np.random.randint(0,len(excluded))
	flip_to_type = excluded[flip_to_ind]
	new_atoms[flip_ind].symbol = flip_to_type
	return new_atoms

def perturb_N_atom(atoms,scale=0.5,fraction=0.25):
	pert_inds = np.random.choice(range(len(atoms)),size=int(len(atoms)*fraction) )
	new_atoms = atoms.copy()
	for pert_ind in pert_inds:
		perturbation = np.random.rand(1,3)[0]
		posneg = 2.*(perturbation - np.min(perturbation))/np.ptp(perturbation)-1
		posneg *= scale
		new_atoms[pert_ind].x += posneg[0]
		new_atoms[pert_ind].y += posneg[1]
		new_atoms[pert_ind].z += posneg[2]
	return new_atoms


def flip_N_atoms(atoms,types,fraction=None):
	fraction = np.random.rand()
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

def mc_step(atoms,types,big_step=False,scale = 0.5,pos_vs_chem_ratio=0.0):
	if big_step:
		step_types = ['perturb_N','flip_N'] 
		step_probs = [pos_vs_chem_ratio, 1-pos_vs_chem_ratio]
	elif not big_step:
		step_types = ['perturb_one','flip_one'] 
		step_probs = [pos_vs_chem_ratio, 1-pos_vs_chem_ratio]
	step_type = np.random.choice(step_types,p=step_probs)
	if step_type == 'flip_N':
		stepped = flip_N_atoms(atoms,types=types)
	elif step_type == 'perturb_N':
		stepped = perturb_N_atom(atoms,scale=scale)
	elif step_type == 'flip_one':
		stepped = flip_one_atom(atoms,types=types)
	elif step_type == 'perturb_one':
		stepped = perturb_one_atom(atoms,scale=scale)
	return stepped

def boltzmann_fac(arr,K=0.5):
	factors = np.exp(-arr/K)
	return factors


target_structure_name = 'supercell_target.cif'
file_prefix = 'iter_%d' % 0
atoms = read(target_structure_name)
write('%s.data' % file_prefix,atoms,format='lammps-data')
full_arr = run_struct(atoms, '%s.data'% file_prefix)
target_arr = full_arr[ : len(atoms), : -1]
target_1 = np.average(target_arr,axis=0)
target_2 = np.var(target_arr,axis=0)

accepted = 1
reduced_structure = read('starting.cif')
file_prefix = 'mc_%d' % 0
write('%s.data' % file_prefix,reduced_structure,format='lammps-data')
full_start_arr = run_struct(reduced_structure, '%s.data'% file_prefix)


Q1,abs_diffs1 = cost_func_1(full_start_arr,target_1,reduced_structure)
Q2,abs_diffs2 = cost_func_2(full_start_arr,target_2,reduced_structure)
print (Q1)
print (Q2)
betas = []


def MC_loop(starting_atoms_in, species, beta, minsteps,moment_weights = (1.0,0.0),n_accept=1,step_size=0.005,target_1=target_1,target_2=target_2,abs_tol = 1.e-14):
	types = [i + 1 for i in range(len(species))]
	count = 1
	hist = []
	full_avg = np.array([])
	last_Q = (1/beta) * 500000
	starting_atoms = starting_atoms_in.copy()
	while count < minsteps and last_Q > abs_tol:
		big_step = np.random.choice([True,False],p=[0.05,0.95])
		try:
			tst_atoms = mc_step(current_atoms,species,big_step=big_step,scale=step_size,pos_vs_chem_ratio=0.0)
			#print (tst_atoms,current_atoms)
		except UnboundLocalError:
			current_atoms = starting_atoms.copy()
			tst_atoms = mc_step(current_atoms,species,big_step=big_step,scale=step_size,pos_vs_chem_ratio=0.0)
		test_comp = get_comp(tst_atoms,species)
		this_prefix = 'temp_%f_iter_%d' % (beta,count)
		write('%s.data' % file_prefix,tst_atoms,format='lammps-data')
		tst_arr = run_struct(tst_atoms, '%s.data'% this_prefix)
		avg_tst = np.average(tst_arr,axis=0)
		if np.shape(full_avg) == (0,):
			full_avg = np.array([avg_tst])
		elif np.shape(full_avg) != (0,):
			full_avg = np.append(full_avg,np.array([avg_tst]),axis=0)
		var_tst = np.average(tst_arr,axis=0)
		
		this_Q1,abs_diffs1 = cost_func_1(tst_arr,target_1,tst_atoms)
		this_Q2,abs_diffs2 = cost_func_2(tst_arr,target_2,tst_atoms)
		test_Q = (moment_weights[0]*this_Q1) + (moment_weights[1]*this_Q2)
		logical_1 = test_Q < last_Q
		rndi = np.random.rand()
		resid_diff = test_Q - last_Q
		bolt1 = np.exp(-resid_diff*beta)
		logical_2 = rndi <= bolt1
		print (count, test_Q, last_Q)
		#print (count,Q1,Q2)
		if logical_1:
			last_Q = test_Q
			current_atoms = tst_atoms
			n_accept += 1
		elif not logical_1 and logical_2:
			last_Q = test_Q
			current_atoms = tst_atoms
			n_accept += 1
		else:
			last_Q = last_Q
			current_atoms = current_atoms
		hist.append(last_Q)
		count += 1
	#print (abs_diffs1)
	#print (abs_diffs2)
	return current_atoms, last_Q, hist, full_avg, n_accept
"""
initial_resid = Q1 + (0.0*Q2)
beta = (1/initial_resid)*10
stepped_atoms, last_Q, resid_hist, n_accept = MC_loop(reduced_structure,['Ca','Mg'],beta,minsteps=100)

print (initial_resid,last_Q)
print (atoms,stepped_atoms)
write('this_step.cif',stepped_atoms)
"""
moment_weights = (1,0.1)
initial_resid = (moment_weights[0]*Q1) + (moment_weights[1]*Q2)
beta_scale = (1/initial_resid)
stepped_atoms=None
#beta_facs=[1.e-3, 1.e1, 1.e3]
beta_facs=[1.e-2,1.e1,1.e2]
for ibeta,beta_fac in enumerate(beta_facs):
	beta = beta_scale* beta_fac
	if stepped_atoms != None:
		reduced_structure = stepped_atoms.copy()
	stepped_atoms, last_Q, resid_hist, full_hist, n_accept = MC_loop(reduced_structure,['Ca','Mg'],beta, minsteps=300, moment_weights=moment_weights)
	np.save('full_avg_hist_%d.npy' % ibeta,np.array(full_hist)) 
	print (initial_resid,last_Q)
	print (atoms,stepped_atoms)
	write('this_step_%d.cif' % ibeta,stepped_atoms)
