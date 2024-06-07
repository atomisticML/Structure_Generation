
from __future__ import print_function
import sys,os
import ctypes
from lib.ga_tools import *
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
    #species = sorted(list(set(list(type_map.values()))))
    species_types = [map_type[spec] for spec in species]
    inp_masses = [masses[typ] for typ in species_types]
    cmds = ["-screen", "none", "-log", "none", "-nocite"]
    lmp = lammps(cmdargs = cmds)

    #def set_atoms(atoms,atid=0):
    #    write('iter_%d.data' % atid,atoms,format='lammps-data')
    #    lmp.command('read_data  iter_%d.data' % atid )
    #    lmp.command('mass  1 180.94788')
    #    lmp.command(f"run {nsteps}")
    import os
    os.path.isfile('%s.data' % fname)
    def run_lammps(dgradflag):

        # simulation settings
        #fname = file_prefix
        lmp.command("clear")
        lmp.command("info all out log")
        lmp.command('units  metal')
        lmp.command('atom_style  atomic')
        #lmp.command("atom_modify    map hash")
        lmp.command('neighbor  2.3 bin')
        # boundary
        #lmp.command('boundary  p p p')
        # read atoms
        lmp.command('read_data  %s.data' % fname )
        for ityp,species_typ in enumerate(species_types):
            imass = inp_masses[ityp]
            lmp.command('mass  %d %f' % (ityp+1,imass))

        # potential settings

        lmp.command(f"pair_style     zero 6.5")
        lmp.command(f"pair_coeff     * *")

        # define compute pace

        if dgradflag:
            lmp.command(f"compute     desc all pace coupling_coefficients.yace 1 1 ")
        else:
            lmp.command(f"compute     desc all pace coupling_coefficients.yace 1 0")

        # run

        lmp.command(f"thermo         1")
        lmp.command(f"thermo_style    custom step temp pe ke etotal c_desc[1][1]")
        #lmp.command(f"run {nsteps}")
        lmp.command(f"timer  off")
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

target_structure_name = 'supercell_target.cif'
file_prefix = 'iter_%d' % 0
atoms = read(target_structure_name)
write('%s.data' % file_prefix,atoms,format='lammps-data')
#full_arr = run_struct(atoms, '%s.data'% file_prefix)
full_arr = run_struct(atoms, file_prefix)
target_arr = full_arr[ : len(atoms), : -1]
target_1 = np.average(target_arr,axis=0)
target_2 = np.var(target_arr,axis=0)

accepted = 1
reduced_structure = read('starting.cif')
file_prefix = 'mc_%d' % 0
write('%s.data' % file_prefix,reduced_structure,format='lammps-data')
full_start_arr = run_struct(reduced_structure, file_prefix)


Q1,abs_diffs1 = cost_func_1(full_start_arr,target_1,reduced_structure)
Q2,abs_diffs2 = cost_func_2(full_start_arr,target_2,reduced_structure)
print (Q1)
print (Q2)
betas = []


def genalg_structure_loop(starting_atoms_in,population_size=4,ngenerations=10):
    countmaxtot = int(population_size*(ngenerations+2))
    q1_weight = 0.8
    q2_weight = 0.2
    igen = 0
    # selection method (only tournament is currently implemented)
    selection_method = 'tournament'

    # cross over (parenting) and mutation hyperparameters
    r_cross = 0.8
    r_mut = 0.2
    # convergence threshold for full function (value of RMSE E + RMSE F at which simulation is terminated" 
    convthr = 0.005
    # fraction of ngenerations to start checking for convergence (convergence checks wont be performed very early)
    conv_check = 1.2
    #best_eval = 999999.999
    conv_flag = False

    all_species = starting_atoms_in.symbols
    cellg = starting_atoms_in.get_cell()
    population = starting_generation(population_size,all_species,cellg,typ='ase')[:-1] + [starting_atoms_in]
    
    def get_cost(tst_atoms,ii,igen):
        this_prefix = 'i_%d_gen_%d' % (ii,igen)
        write('%s.data' % this_prefix,tst_atoms,format='lammps-data')
        write('%s.cif' % this_prefix,tst_atoms,format='cif')

        tst_arr = run_struct(tst_atoms, this_prefix)
        this_Q1,abs_diffs1 = cost_func_1(tst_arr,target_1,tst_atoms)
        this_Q2,abs_diffs2 = cost_func_2(tst_arr,target_2,tst_atoms)
        return (q1_weight*this_Q1)+(q2_weight*this_Q2)

    scores = [get_cost(test_atoms,ii,igen) for ii,test_atoms in enumerate(population)]

    best_eval = min(scores)
    best_id = scores.index(best_eval)

    best_evals = [best_eval]
    best_inds = ['i_%d_gen_%d' % (best_id,igen)]
    while igen <= ngenerations and best_eval > convthr and not conv_flag:
        slct = Selector(selection_style = selection_method)
        selected = [slct.selection(population, scores) for creature_idx in range(population_size)]
        del slct
        children = list()
        for ii in range(0,population_size,2):
            p1,p2 = selected[ii], selected[ii+1]
            # crossover and mutation
            rndcross, rndmut = tuple(np.random.rand(2).tolist())
            if rndcross <= r_cross:
                cs = crossover(p1,p2)
            else:
                cs = [p1,p2]
            for c in cs:
                # mutation
                if rndmut <= r_mut:
                    muttyp = mutation_type_from_prob(choice_probs = {'perturb_one': 0.5, 'perturb_N' : 0.5, 'flip_one' : 0.5, 'flip_N' : 0.5})
                    c = mutation(c,muttyp,types=['Ca','Mg'],scale=1.5)
                children.append(c)
        igen += 1
        population = children
        scores = [get_cost(test_atoms,ii,igen) for ii,test_atoms in enumerate(population)]
        #print ('gen score',scores)
        print('best score for gen %d: %f' % (igen,min(scores)))

        best_eval = min(scores)
        best_id = scores.index(best_eval)
        best_evals.append(best_eval)
        best_inds.append('i_%d_gen_%d' % (best_id,igen) )
    print (best_evals)
    print (best_inds)

genalg_loop(reduced_structure,population_size=60,ngenerations=100)
