import numpy as np
from ase.io import read,write

def generate_random_integers(sum_value, n):
    # Generate a list of n - 1 random integers
    random_integers = [np.random.randint(0, sum_value) for _ in range(n - 1)]
    # Ensure the sum of the random integers is less than or equal to sum_value
    random_integers.sort()
    # Calculate the Nth integer to ensure the sum is equal to sum_value
    random_integers.append(sum_value - sum(random_integers))
    return random_integers

def starting_generation(pop_size,all_species,cell,typ='ase',nchem = 1):
    pop = []
    if typ == 'ase':
        from ase.io import read,write
        from ase import Atoms,Atom
        from ase.ga.utilities import closest_distances_generator, CellBounds
        from ase.ga.startgenerator import StartGenerator
        from ase.data import atomic_numbers
        volume = np.dot(cell[2],np.cross(cell[0],cell[1]))
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
        return pop

def at_to_lmp(atoms,index):
    #s=generate.format(xx,yy,zz,dx,dz,dz,n_atoms,seed,index,random.uniform(0.0,1.0))
    fname = 'ats_%s.data' % index
    write(fname,atoms,format='lammps-data')
    generate=\
"""
units           metal
boundary        p p p

read_data  {}

log log_{}.lammps
mass 1 58.6934

pair_style hybrid/overlay soft 2.0 mliap model mliappy LATER descriptor ace coupling_coefficients.yace
pair_coeff * * soft 100
pair_coeff * * mliap Ni

thermo 10
fix nve all nve
fix lan all langevin 5000 100 1.0 48279

velocity all create 10000 4928459 dist gaussian
"""
    #s = generate.format(fname,index,np.random.uniform(0.0,1.0))
    s = generate.format(fname,index)
    return s

#pair_style hybrid/scaled {} soft 2.0 1.0 LATER mliap model mliappy LATER descriptor pace coupling_coefficients.yace 
#def internal_generate_cell(atoms,index):
def internal_generate_cell(index):
    from ase.build import bulk
    template = bulk('Ni',cubic=True)*(1,2,2)
    all_species = template.symbols
    cellg = template.get_cell()
    rnd = starting_generation(1,all_species,cellg,typ='ase')[0]
    s = at_to_lmp(rnd,index)
    return s

#internal_generate_cell(0)
