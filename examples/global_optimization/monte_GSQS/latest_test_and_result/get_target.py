
from __future__ import print_function
import sys,os
import ctypes
import numpy as np
from ase.io import read,write
from lammps import lammps, LMP_TYPE_ARRAY, LMP_STYLE_GLOBAL

# get MPI settings from LAMMPS
def run_struct(atoms,fname):
		
	lmp = lammps()
	me = lmp.extract_setting("world_rank")
	nprocs = lmp.extract_setting("world_size")


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
		lmp.command("boundary	p p p")
		lmp.command("atom_modify	map hash")
		lmp.command('neighbor  2.3 bin')
		# boundary
		lmp.command('boundary  p p p')
		# read atoms
		lmp.command('read_data  %s.data' % fname )
		lmp.command('mass  1 40.078')
		lmp.command('mass  2 24.305')

		# potential settings

		lmp.command(f"pair_style 	zero 4.5")
		lmp.command(f"pair_coeff 	* *")

		# define compute pace

		if dgradflag:
			lmp.command(f"compute 	pace all pace coupling_coefficients.yace 1 1")
		else:
			lmp.command(f"compute 	pace all pace coupling_coefficients.yace 1 0")

		# run

		lmp.command(f"thermo 		100")
		lmp.command(f"run {nsteps}")

	# declare simulation/structure variables

	nsteps = 0
	ntypes = 1

	# declare compute pace variables

	bikflag = 1

	# NUMBER of descriptors
	#nd = 236 
	#nd = 73

	# run lammps with dgradflag on


	dgradflag = 0

	run_lammps(dgradflag)

	# get global snap array
	lmp_pace = lmp.numpy.extract_compute("pace", LMP_STYLE_GLOBAL, LMP_TYPE_ARRAY)
	#print (lmp_pace)
	#print (np.shape(lmp_pace))
	#for row in lmp_pace:
	#	print (row)
	descriptor_grads = lmp_pace[ : len(atoms), : -1]


	#del lmp
	return descriptor_grads

import sys

#f = sys.argv[-1]

import glob

start = 'supercell_target.cif'
file_prefix = 'iter_%d' % 0
atoms = read(start)

write('%s.data' % file_prefix,atoms,format='lammps-data')
start_arr = run_struct(atoms, '%s.data'% file_prefix)
print (start_arr)
print (np.shape(start_arr))
avg_start = np.average(start_arr,axis=0)
print (np.shape(avg_start))

var_start = np.var(start_arr,axis=0)
np.save('target_descriptors.npy',avg_start)
np.save('target_var_descriptors.npy',var_start)
