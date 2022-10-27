"""
compute_pace.py
"""

from __future__ import print_function
import sys, os
import ctypes
import numpy as np
from ase.io import read,write
from lammps import lammps, LMP_TYPE_ARRAY, LMP_STYLE_GLOBAL

# get MPI settings from LAMMPS
def run_struct(f):

	# runs and saves 1st and 2nd moments of ace descriptors within a structure

	# set to True to save the full, per atom descriptor array
	save_full_array = False

	# set to false if you do not want to calculate skewness or kurtosis from scipy

	sci_vals = True

	#dictionary of masses (example for H, O)
	mass_per_lammps_type = { 
	1: 1.00,
	2: 15.999
	}

	file_prefix = f.split('.')[0]

	atoms = read(f)
	natoms = len(atoms)
	lmp = lammps()

	me = lmp.extract_setting("world_rank")
	nprocs = lmp.extract_setting("world_size")

	write('%s.data' % file_prefix,atoms,format='lammps-data')

	cmds = ["-screen", "none", "-log", "none"]
	lmp = lammps(cmdargs = cmds)

	def run_lammps():

		# simulation settings
		fname = file_prefix
		lmp.command("clear")
		lmp.command("info all out log")
		lmp.command('units  metal')
		lmp.command('atom_style  atomic')
		lmp.command("boundary	p p p")
		lmp.command("atom_modify	map hash")
		lmp.command('neighbor  2.3 bin')
		# boundary
		lmp.command('boundary  p p p')
		# read atoms
		lmp.command('read_data  %s.data' % fname )
		# define lammps types and corresponding masses
		for typ in mass_per_lammps_type.keys():
			lmp.command('mass  %d %f' % (typ,mass_per_lammps_type[typ]))


		lmp.command(f"pair_style 	zero 4.55")
		lmp.command(f"pair_coeff 	* *")
		# define compute pace

		lmp.command(f"compute 	pace all pace coupling_coefficients.yace 1 0")
		# run
		lmp.command(f"thermo 		1")
		lmp.command(f"run 0")

	run_lammps()

	# get global pace array
	lmp_pace = lmp.numpy.extract_compute("pace", LMP_STYLE_GLOBAL, LMP_TYPE_ARRAY)
	# extract the first 'natoms' rows, corresponding to descriptor values for each atom i=0,1 ... natoms (exclusive)
	eis = lmp_pace[ : natoms, : -1]
	if save_full_array:
		np.save('%s_descs.npy' % file_prefix,eis)
	elif not save_full_array:
		np.save('%s_avg.npy' % file_prefix, np.average(eis,axis=0))
		np.save('%s_var.npy' % file_prefix, np.var(eis,axis=0))
		if sci_vals:
			from scipy import stats
			np.save('%s_skew.npy' % file_prefix, stats.skew(eis,axis=0))
			np.save('%s_kurt.npy' % file_prefix, stats.kurtosis(eis,axis=0))
	lmp.close()
	return None

import glob

#example calculates descriptors for all structure files with the 'xyz' suffix
# It will read files through ASE, and any ASE-compatible structure file may be used
for f in sorted(glob.glob('*.xyz')):
	print ('running %s' % f)
	run_struct(f)
