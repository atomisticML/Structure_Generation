import numpy as vnp
from ase.io import read,write

#manually tabulated minimum soft core cutoffs (if not tabulated here ASE-generated default will be used)
min_soft = {'Cr':1.65,'Fe':1.7,'Si':1.2,'V':1.65}
# TODO - make minimum cutoffs an input


def build_target(start):
    from ase.io import read,write
    from ase import Atoms,Atom
    from ase.ga.utilities import closest_distances_generator, CellBounds
    from ase.ga.startgenerator import StartGenerator
    from ase.data import atomic_numbers, atomic_names, atomic_masses, covalent_radii
    #from __future__ import print_function
    import sys,os
    import ctypes
    import numpy as np
    from lammps import lammps, LMP_TYPE_ARRAY, LMP_STYLE_GLOBAL

    # get mpi settings from lammps
    def run_struct(atoms,fname,maxcut=6.0):
            
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
            utypes = []
            for atom in atoms:
                if atom.symbol not in utypes:
                    utypes.append(atom.symbol)
            for ind,typ in enumerate(utypes):
                number = atomic_numbers[typ]
                mass = atomic_masses[number]
                lmp.command('mass   %d %f' % (ind+1,mass))

            lmp.command("pair_style 	zero %f" % maxcut)
            lmp.command(f"pair_coeff 	* *")

            if dgradflag:
                lmp.command(f"compute 	pace all pace coupling_coefficients.yace 1 1")
            else:
                lmp.command(f"compute 	pace all pace coupling_coefficients.yace 1 0")

            # run

            lmp.command(f"thermo 		100")
            #lmp.command(f"run {nsteps}")
            lmp.command(f"run 0")


        # declare compute pace variables

        dgradflag = 0
        run_lammps(dgradflag)
        lmp_pace = lmp.numpy.extract_compute("pace", LMP_STYLE_GLOBAL, LMP_TYPE_ARRAY)
        descriptor_grads = lmp_pace[ : len(atoms), : -1]
        return descriptor_grads
    #start = 'supercell_target.cif'
    file_prefix = 'iter_%d' % 0
    atoms = read(start)

    write('%s.data' % file_prefix,atoms,format='lammps-data')
    start_arr = run_struct(atoms, '%s.data'% file_prefix)
    avg_start = vnp.average(start_arr,axis=0)
    var_start = vnp.var(start_arr,axis=0)
    vnp.save('target_descriptors.npy',avg_start)
    vnp.save('target_var_descriptors.npy',var_start)
    return avg_start, var_start

def rand_comp(model_dict):
    old_vals = vnp.array(list(model_dict.values()))
    vals = vnp.random.rand(len(old_vals))
    newvals = vals/vals.sum()
    new_dict = {list(model_dict.keys())[ii]:newvals[ii] for ii in range(len(newvals))}
    return new_dict
    
def generate_random_integers(sum_value, n):
    # Generate a list of n - 1 random integers
    random_integers = [vnp.random.randint(0, sum_value) for _ in range(n - 1)]
    # Ensure the sum of the random integers is less than or equal to sum_value
    random_integers.sort()
    # Calculate the Nth integer to ensure the sum is equal to sum_value
    random_integers.append(sum_value - sum(random_integers))
    return random_integers

def starting_generation(pop_size,all_species,cell,typ='ase',nchem = 1,use_template=None):
    pop = []
    if typ == 'ase':
        from ase.io import read,write
        from ase import Atoms,Atom
        from ase.ga.utilities import closest_distances_generator, CellBounds
        from ase.ga.startgenerator import StartGenerator
        from ase.data import atomic_numbers
        volume = vnp.dot(cell[2],vnp.cross(cell[0],cell[1]))
        # Target cell volume for the initial structures, in angstrom^3
        #volume = 240.        # Specify the 'building blocks' from which the initial structures
        # will be constructed. Here we take single Ag atoms as building
        # blocks, 24 in total.
        #blocks = [('Ag', 24)]
        ## We may also write:
        #blocks = ['Ag'] * 24
        dsize = int(pop_size/nchem)
        sort_specs = sorted(all_species)
        nats = len(all_species)
        uspecs = list(set(sort_specs))
        if len(uspecs) == 1:
            block_sets = [ [(uspecs[0],len(sort_specs))] ]* nchem
        else:
            block_sets = []
            for iii in range(nchem):
                natoms_per_type = generate_random_integers(len(sort_specs), len(uspecs))
                block_sets.append( [ tuple([uspecs[ij],natoms_per_type[ij]]) for ij in range(len(uspecs)) ] )
            #blocks = [('Ti', 4), ('O', 8)] 
        # Generate a dictionary with the closest allowed interatomic distances
        Zs = [ atomic_numbers[sym] for sym in list(set(all_species))]
        blmin = closest_distances_generator(atom_numbers=Zs, ratio_of_covalent_radii=0.5)
        natoms = len(all_species)        
        slab = Atoms('',cell=cell, pbc=True)
        for block in block_sets:
            # Initialize the random structure generator
            #sg = StartGenerator(slab, all_species, blmin, box_volume=volume,
            sg = StartGenerator(slab, block, blmin, number_of_variable_cell_vectors=0) 
            # and add them to the database
            for i in range(dsize):
                a = sg.get_new_candidate()
                pop.append(a)
            if use_template:
                pop = [use_template] + pop
        return pop

def at_to_lmp(atoms,index,temperature=10000.0,min_typ='temp', coefftypes=True,soft_strength=10000.0):
    from ase.data import atomic_numbers, atomic_names, atomic_masses, covalent_radii
    from ase.ga.utilities import closest_distances_generator
    #s=generate.format(xx,yy,zz,dx,dz,dz,n_atoms,seed,index,random.uniform(0.0,1.0))
    fname = 'ats_%s.data' % index
    write(fname,atoms,format='lammps-data')
    generate=\
"""
units           metal
boundary        p p p

read_data  {}

log log_{}.lammps
""" 
    s1 = generate.format(fname,index)
    atutypes =[]
    for atom in atoms:
        if atom.symbol not in atutypes:
            atutypes.append(atom.symbol)
    if coefftypes:
        utypes = get_desc_count(coefffile='coupling_coefficients.yace',return_elems=True)
    else:
        utypes = atutypes.copy()
    if len(atutypes) < len(utypes):
        utypes = atutypes.copy()
    typstr = ' ' 
    try:
        Zs = [ atomic_numbers[sym] for sym in utypes]
        blmin = closest_distances_generator(atom_numbers=Zs, ratio_of_covalent_radii=0.8)
        min_soft_utypes = list(blmin.values())
    except KeyError:
        print('missing tabulated minimum bonds, applying minimums from ASE')
    for indtyp,utype in enumerate(utypes):
        atnum = atomic_numbers[utype]
        atmass = atomic_masses[atnum]
        s1 += 'mass %d %f\n' % (indtyp+1, atmass)
        typstr += ' %s' % utype
    #print('typstr final: %s' % typstr)
    #print (utypes,atutypes)
    if min_typ =='temp':
        generate2a=\
"""pair_style hybrid/overlay soft %2.3f mliap model mliappy LATER descriptor ace coupling_coefficients.yace
pair_coeff * * soft %f\n""" % (max(min_soft_utypes),soft_strength)
        generate2b=\
"""pair_coeff * * mliap %s

thermo 10
fix nve all nve
fix lan all langevin %d 100 1.0 48279

velocity all create %d 4928459 dist gaussian
""" % (typstr, int(temperature), int(temperature*2))
        generate2 = generate2a+generate2b
    elif min_typ == 'box':
        generate2a =\
"""pair_style hybrid/overlay soft %2.3f mliap model mliappy LATER descriptor ace coupling_coefficients.yace
pair_coeff * * soft %f\n""" % (max(min_soft_utypes),soft_strength)
        generate2b =\
"""pair_coeff * * mliap %s

thermo 10
fix 1b all box/relax iso 0.0 vmax 0.001
velocity all create 0.0001 4928459 dist gaussian""" % typstr
        generate2 = generate2a+generate2b
    else:
        generate2a =\
"""pair_style hybrid/overlay soft %2.3f mliap model mliappy LATER descriptor ace coupling_coefficients.yace
pair_coeff * * soft %f\n""" % (max(min_soft_utypes),soft_strength)
        generate2b = \
"""pair_coeff * * mliap %s

thermo 10
velocity all create 0.0001 4928459 dist gaussian""" % typstr
        generate2 = generate2a+generate2b
#minimize 1e-8 1e-8 1000 1000"
    #s = generate.format(fname,index,vnp.random.uniform(0.0,1.0))
    s = s1 + generate2
    return s

from hnf import *
def internal_generate_cell(index,desired_size=4,template=None,desired_comps={'Ni':1.0},use_template=None,min_typ='temp',soft_strength=10000):
    if template == None:
        from ase.build import bulk
        from ase import Atoms,Atom
        chems = list(desired_comps.keys())
        template = Atoms([chems[0]]*desired_size)
        atoms_base = bulk(chems[0])
        vol_base = vnp.dot(vnp.cross(atoms_base.get_cell()[0],atoms_base.get_cell()[1]),atoms_base.get_cell()[2])
        a_simp = vol_base**(1/3)
        #cells = get_hnfs(hnf_trs=[desired_size])
        cells_all = get_hnfs(hnf_trs=[desired_size])
        toli = int((desired_size)**(1/3))
        cells= limit_mats_len(cells_all,desired_size,tol=0.16)
        #print ('ncells',len(cells))
        try:
            cell = a_simp * cells[vnp.random.choice(range(len(cells)))]
        except ValueError:
            cell = a_simp * cells_all [vnp.random.choice(range(len(cells_all)))]
        #print ('hnf',cell)
        norms = np.array([np.linalg.norm(k) for k in cell])
        norms /= (desired_size)
        #print('hnf norms',norms)
        #cell = a_simp * cells[-2]
        template.set_cell(cell)
        template.set_scaled_positions(vnp.random.uniform(0,1,(desired_size,3)))
    else:
        tempalte = template
    new_comps = {elem:int(round(len(template)*cmp))/len(template) for elem,cmp in desired_comps.items()}
    print ('for structure of size:%d'% len(template),'desired compositions:', desired_comps,'will be replaced with', new_comps)
    all_species = get_target_comp_s(desired_size=len(template),desired_comps=new_comps,parent_cell=template.get_cell())
    #all_species = get_target_comp(desired_size=len(template),desired_comps=new_comps,parent_cell=template.get_cell())
    if len(all_species) >= len(template):
        print('doing cheap fix for cell size')
        all_species=all_species[:len(template)]
    elif len(all_species) <= len(template):
        print('doing cheap grow for cell size')
        diff = len(template)-len(all_species)
        all_species = all_species + all_species[:diff]
    #print ('all specs vs template',len(all_species),len(template))
    assert len(all_species)== len(template), "composition list size must match size of template atoms"
    cellg = template.get_cell()
    if use_template:
        rnd = starting_generation(1,all_species,cellg,typ='ase',use_template=template)[0]
    else:
        rnd = starting_generation(1,all_species,cellg,typ='ase')[0]
    s = at_to_lmp(rnd,index,min_typ=min_typ,soft_strength=soft_strength)
    return s

def get_comp(atoms,symbols):
    comps = {symbol: 0.0 for symbol in symbols}
    counts = {symbol: 0 for symbol in symbols}
    atsymbols = [atom.symbol for atom in atoms]
    for atsymbol in atsymbols:
        counts[atsymbol] +=1
    for symbol in symbols:
        comps[symbol] = counts[symbol]/len(atoms)
    return comps

def flip_one_atom(atoms,types):
    new_atoms = atoms.copy()
    flip_ind = vnp.random.randint(0,len(atoms))
    flip_current = new_atoms[flip_ind].symbol
    excluded = [typ for typ in types if typ != flip_current]
    flip_to_ind = vnp.random.randint(0,len(excluded))
    flip_to_type = excluded[flip_to_ind]
    new_atoms[flip_ind].symbol = flip_to_type
    return new_atoms

def flip_N_atoms(atoms,types,fraction=0.25):
    pert_inds = vnp.random.choice(range(len(atoms)),size=int(len(atoms)*fraction) )
    new_atoms = atoms.copy()
    for pert_ind in pert_inds:
        flip_ind = vnp.random.randint(0,len(atoms))
        flip_current = new_atoms[flip_ind].symbol
        excluded = [typ for typ in types if typ != flip_current]
        flip_to_ind = vnp.random.randint(0,len(excluded))
        flip_to_type = excluded[flip_to_ind]
        new_atoms[flip_ind].symbol = flip_to_type
    return new_atoms


def add_atom(atoms,symbols,tol = 0.5):
    from ase import Atom,Atoms
    from ase.ga.utilities import closest_distances_generator
    from ase.data import atomic_numbers
    from ase.neighborlist import primitive_neighbor_list
    blmin = closest_distances_generator(atom_numbers=[atomic_numbers[symbol] for symbol in symbols] + [atomic_numbers['Ne']], ratio_of_covalent_radii=0.5)
    def readd():
        symbol = vnp.random.choice(symbols)
        rnd_pos_scale = vnp.random.rand(1,3)
        rnd_pos = vnp.matmul(atoms.get_cell(),rnd_pos_scale.T)
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
    conds = all([ vnp.linalg.norm(tst_dist) >=  blmin[(atomic_numbers[symbol] , tst_bonds[i][1])] for i,tst_dist in enumerate(tst_dists)])
    while not conds:
        symbol , tst_atoms, at_dists , rnd_pos, bond_types = readd()
        syms = [tst_atom.symbol for tst_atom in tst_atoms]
        tst_id = syms.index('Ne')
        tst_dists = at_dists[tst_id]
        tst_bonds = bond_types[tst_id]
        #conds = all([ vnp.linalg.norm(tst_dist) >= tol for tst_dist in tst_dists])
        #conds = all([ vnp.linalg.norm(tst_dist) >= blmin[tst_bonds[i]] for i,tst_dist in enumerate(tst_dists)])
        conds = all([ vnp.linalg.norm(tst_dist) >=  blmin[(atomic_numbers[symbol] , tst_bonds[i][1])]-tol for i,tst_dist in enumerate(tst_dists)])
    atoms.append(Atom(symbol,rnd_pos))
    return atoms

def get_target_comp_s(desired_size,desired_comps,parent_cell):
    num_ats = {key:int(desired_comps[key]*desired_size) for key in list(desired_comps.keys())}
    symbols =[]
    for key,nrepeat in num_ats.items():
        symbols.extend([key]*nrepeat)
    return symbols

def get_target_comp(desired_size,desired_comps,parent_cell):
    from ase import Atoms
    initial_a = Atoms('',pbc=True)
    initial_a.set_cell(parent_cell)

    current_comps = {key:0.0 for key in list(desired_comps.keys())}
    current_size = 0

    chems = list(desired_comps.keys())
    comp_conds = all([current_comps[chem] == desired_comps[chem] for chem in chems])
    these_atoms = initial_a.copy()
    max_iter = 100
    this_iter = 0 
    toli = 1/(desired_size)
    while current_size <= desired_size or not comp_conds and this_iter <= max_iter:
        #t = 1000 * time.time() # current time in milliseconds
        #vnp.random.seed(int(t) % 2**32)
        if current_size < desired_size:
            tst_ats = add_atom(these_atoms,chems)
            tst_comps =  get_comp(tst_ats, chems)
            comp_conds = all([tst_comps[chem] == desired_comps[chem] for chem in chems])
            these_atoms = tst_ats.copy()
            current_size = len(these_atoms)
            Qi = [vnp.abs(tst_comps[chem] - desired_comps[chem]) for chem in chems]
            Qi = round(vnp.sum(Qi),8)
        elif current_size == desired_size:
            #tst_ats = flip_one_atom(these_atoms,chems)
            tst_ats = flip_N_atoms(these_atoms,chems,fraction=vnp.random.rand())
            tst_comps =  get_comp(tst_ats, chems)
            comp_conds = all([tst_comps[chem] == desired_comps[chem] for chem in chems])
            tstQi = [vnp.abs(tst_comps[chem] - desired_comps[chem]) for chem in chems]
            tstQi = round(vnp.sum(tstQi),8)
            #print (len(these_atoms),tst_comps,tstQi)
            if tstQi <= Qi+toli:
                these_atoms = tst_ats.copy()
                Qi = tstQi
                #print ('in composition loop',tst_comps,comp_conds)
        if current_size==desired_size and comp_conds:
            break
        this_iter += 1
    return list(tst_ats.symbols)

def get_desc_count(coefffile,return_elems=False):
    with open(coefffile,'r') as readcoeff:
        lines = readcoeff.readlines()
        elemline = [line for line in lines if 'elements' in line][0]
        elemstr = elemline.split(':')[-1]
        elemstr2 = elemstr.replace('[','')
        elemstr3 = elemstr2.replace(']','')
        elemstr4 = elemstr3.replace(',','')
        elems = elemstr4.split()
        nelements = len(elems)
        desclines = [line for line in lines if 'mu0' in line]
    if not return_elems:
        return int(len(desclines)/nelements)
        #return int(len(desclines))
    else:
        return elems

#rnd = internal_generate_cell(0,desired_comps={'W':0.9,'H':0.1})
#rnd = internal_generate_cell(index,desired_comps={'W':0.9,'H':0.1})
#rnd = internal_generate_cell(1,desired_comps={'Cr':1.0})
#print(rnd)
