import numpy as np
#from ase.io import read,write
#from ase.ga.utilities import closest_distances_generator, CellBounds
#from ase.ga.startgenerator import StartGenerator

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
        #volume = 240.

        # Specify the 'building blocks' from which the initial structures
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
        blmin = closest_distances_generator(atom_numbers=Zs,
                                            ratio_of_covalent_radii=0.5)

        natoms = len(all_species)

        slab = Atoms('',cell=cell, pbc=True)

        for block in block_sets:
            # Initialize the random structure generator
            #sg = StartGenerator(slab, all_species, blmin, box_volume=volume,
            sg = StartGenerator(slab, block, blmin,
                                number_of_variable_cell_vectors=0) #splits=

            # Generate N random structures
            # and add them to the database
            for i in range(dsize):
                a = sg.get_new_candidate()
                pop.append(a)
        return pop

#pop_size = 20
#cell = np.eye(3)*6.3
#all_species = ['Ca']*12 + ['Mg']*12
#p = starting_generation(pop_size,all_species,cell,typ='ase')
#print (p)

# tournament selection
def tournament_selection(population, scores, k=3, inputseed=None):
    if inputseed != None:
        np.random.seed(inputseed)
    selection_ix = np.random.randint(len(population))
    for ix in np.random.randint(0, len(population), k-1):
        # check if better (e.g. perform a tournament)
        if scores[ix] < scores[selection_ix]:
            selection_ix = ix
    return population[selection_ix]


class Selector:
    def __init__(self,selection_style = 'tournament'):
        self.selection_style = selection_style
        self.set_selector()

    def set_selector(self):
        if self.selection_style == 'tournament':
            self.selection = tournament_selection

# simple crossover function for atomic coordinates
def crossover(p1, p2, inputseed = None):
    assert len(p1) == len(p2), "parents must have the same length"
    psize = len(p1)
    if inputseed != None:
        np.random.seed(inputseed)
    cross_point = np.random.randint(1, psize-1)
    c1 = p1[:cross_point] + p2[cross_point:]
    c2 = p2[:cross_point] + p1[cross_point:]

    return [c1, c2]


def perturb_one_atom(atoms,scale=0.5,max_attempt=100,apply_to='ase'):
    if apply_to == 'ase':
        from ase.ga.utilities import closest_distances_generator
        from ase.data import atomic_numbers
        import ase
        import ase.neighborlist
        cutoffs = ase.neighborlist.natural_cutoffs(atoms)
        sym_num_map = {sym:atomic_numbers[sym] for sym in atoms.symbols}
        #nl = ase.neighborlist.neighbor_list('ijd', atoms, max(cutoffs))
        nl = ase.neighborlist.neighbor_list('ijd', atoms, cutoffs)
        #nl.update(atoms)
        #indices, offsets = nl.get_neighbors(0)`
        Zs = [ atomic_numbers[sym] for sym in list(set(atoms.symbols))]
        blmin = closest_distances_generator(atom_numbers=Zs,
                                            ratio_of_covalent_radii=0.5)
        good_pert = False
        nattempt =0
        new_atoms = atoms.copy()
        while not good_pert and nattempt < max_attempt:
            new_atoms = atoms.copy()
            pert_ind = np.random.randint(0,len(atoms))
            perturbation = np.random.rand(1,3)[0]
            posneg = 2.*(perturbation - np.min(perturbation))/np.ptp(perturbation)-1
            posneg *= scale
            new_atoms[pert_ind].x += posneg[0]
            new_atoms[pert_ind].y += posneg[1]
            new_atoms[pert_ind].z += posneg[2]
            nl = ase.neighborlist.neighbor_list('ijd', new_atoms, cutoffs)
            sym_num_map = {sym:atomic_numbers[sym] for sym in new_atoms.symbols}
            pair_dist_flags = [ blmin[tuple(sorted([sym_num_map[new_atoms.symbols[i]], sym_num_map[ new_atoms.symbols[nl[1][i_ind]]] ]))] > nl[2][i_ind] for i_ind, i in enumerate(nl[0])]
            if not any(pair_dist_flags):
                good_pert = True
            nattempt += 1
        return new_atoms

        """
        new_atoms = atoms.copy()
        pert_ind = np.random.randint(0,len(atoms))
        perturbation = np.random.rand(1,3)[0]
        posneg = 2.*(perturbation - np.min(perturbation))/np.ptp(perturbation)-1
        posneg *= scale
        new_atoms[pert_ind].x += posneg[0]
        new_atoms[pert_ind].y += posneg[1]
        new_atoms[pert_ind].z += posneg[2]
        return new_atoms
        """

    elif apply_to == 'raw_positions':
        new_atoms = atoms.copy()
        pert_ind = np.random.randint(0,len(atoms))
        perturbation = np.random.rand(1,3)[0]
        posneg = 2.*(perturbation - np.min(perturbation))/np.ptp(perturbation)-1
        posneg *= scale
        new_atoms[pert_ind][0] += posneg[0]
        new_atoms[pert_ind][1] += posneg[1]
        new_atoms[pert_ind][2] += posneg[2]
        return new_atoms

#from ase import Atoms,Atom
#atoms = Atoms(['Ag']*20,positions=np.random.uniform(-8,8, (20,3)))
#atoms.set_cell(np.eye(3)*16)
#perturb_one_atom(atoms)


def flip_one_atom(atoms,types):
    new_atoms = atoms.copy()
    flip_ind = np.random.randint(0,len(atoms))
    flip_current = new_atoms[flip_ind].symbol
    excluded = [typ for typ in types if typ != flip_current]
    flip_to_ind = np.random.randint(0,len(excluded))
    flip_to_type = excluded[flip_to_ind]
    new_atoms[flip_ind].symbol = flip_to_type
    return new_atoms


def perturb_N_atoms(atoms,scale=0.5,max_attempt = 100, fraction=0.25):
    from ase.ga.utilities import closest_distances_generator
    from ase.data import atomic_numbers
    import ase
    import ase.neighborlist
    cutoffs = ase.neighborlist.natural_cutoffs(atoms)
    sym_num_map = {sym:atomic_numbers[sym] for sym in atoms.symbols}
    #nl = ase.neighborlist.neighbor_list('ijd', atoms, max(cutoffs))
    nl = ase.neighborlist.neighbor_list('ijd', atoms, cutoffs)
    #nl.update(atoms)
    #indices, offsets = nl.get_neighbors(0)`
    Zs = [ atomic_numbers[sym] for sym in list(set(atoms.symbols))]
    blmin = closest_distances_generator(atom_numbers=Zs,
                                        ratio_of_covalent_radii=0.5)

    good_pert = False
    nattempt =0
    new_atoms = atoms.copy()
    while not good_pert and nattempt < max_attempt:
        pert_inds = np.random.choice(range(len(atoms)),size=int(len(atoms)*fraction) )
        for pert_ind in pert_inds:
            new_atoms = atoms.copy()
            perturbation = np.random.rand(1,3)[0]
            posneg = 2.*(perturbation - np.min(perturbation))/np.ptp(perturbation)-1
            posneg *= scale
            new_atoms[pert_ind].x += posneg[0]
            new_atoms[pert_ind].y += posneg[1]
            new_atoms[pert_ind].z += posneg[2]
        nl = ase.neighborlist.neighbor_list('ijd', new_atoms, cutoffs)
        sym_num_map = {sym:atomic_numbers[sym] for sym in new_atoms.symbols}
        pair_dist_flags = [ blmin[tuple(sorted([sym_num_map[new_atoms.symbols[i]], sym_num_map[ new_atoms.symbols[nl[1][i_ind]]] ]))] > nl[2][i_ind] for i_ind, i in enumerate(nl[0])]
        if not any(pair_dist_flags):
            good_pert = True
        nattempt += 1
    if not good_pert:
        print ("WARNING: this mutation has a bad neighbor distance - try changing your scale parameter")
    return new_atoms


#from ase import Atoms,Atom
#atoms = Atoms(['Ag']*20,positions=np.random.uniform(-8,8, (20,3)))
#atoms.set_cell(np.eye(3)*16)
#perturb_N_atoms(atoms)


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


def mutation_type_from_prob(choice_probs = {'perturb_one': 0.75, 'perturb_N' : 0.25, 'flip_one' : 0.0, 'flip_N' : 0.0}):
    choices = list(choice_probs.keys())
    probs = list(choice_probs.values())
    mut_typ = np.random.choice(choices,p=probs)
    return mut_typ


def mutation(current_atoms,mutation_type='perturb_N',types=['Ag'],scale=0.5):
    mutation_types = {
    'perturb_one' : perturb_one_atom,
    'perturb_N' : perturb_N_atoms,
    'flip_one' : flip_one_atom,
    'flip_N' : flip_N_atoms,
    }
    if 'flip' in mutation_type:
        return mutation_types[mutation_type](current_atoms,types)
    else:
        return mutation_types[mutation_type](current_atoms,scale)
