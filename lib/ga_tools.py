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

def generate_occs(target, natoms):
    if type(target) == dict:
        counts = {k: round(v * natoms) for k, v in target.items()}
        diff = natoms - sum(counts.values())
        while diff != 0:
            for k in counts:
                if diff > 0:
                    counts[k] += 1
                    diff -= 1
                elif diff < 0 and counts[k] > 0:
                    counts[k] -= 1
                    diff += 1
        occs = [k for k, v in counts.items() for _ in range(v)]
        np.random.shuffle(occs)
    elif type(target) == list:
        occs = np.random.choice(target,natoms).tolist()
    return occs

def starting_generation(pop_size,all_species,cell,typ='ase',nchem = 1,**kwargs):
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

    elif typ == 'lattice':
        from ase.build import bulk
        try:
            parent = kwargs['parent']
            sc = kwargs['s'] # ( 2,2,1 ) or some other tuple that defines supercell
            #atoms = parent*sc
            occ_concs = kwargs['occupation_concentrations'] # {'Mg':0.25,'Ca':0.75} or ['Ca','Mg'] for random sampling
        except KeyError:
            try:
                typ = kwargs['type'] # a SINGLE parent element: e.g. Mg OR Ta OR Ni ...
                lattice = kwargs['lattice'] # 'bcc','fcc',... (just cubic for now to keep it simple)
                a = kwargs['a'] # lattice constant
                #syms = [s for s in atoms.symbols]
                sc = kwargs['s'] # ( 2,2,1 ) or some other tuple that defines supercell
                occ_concs = kwargs['occupation_concentrations']
            except KeyError:
                raise ValueError(" need to supply the type, lattice, lattice constant, and supercell factor to use 'lattice' generation type: example kwargs={ 'type':'Mg', 'lattice':'bcc', 'a':4.25, 's':(2,2,2)}")
            assert lattice in ['bcc','fcc','sc'], "must have cubic lattice in starting_generation=lattice type"
            parent = bulk(typ,lattice,a,cubic=True)
            atoms = parent*sc
            for i in range(pop_size):
                occs = generate_occs(occ_concs, len(atoms))
                print ('occs in initial pop %d' %i, occs)
                newatoms = atoms.copy()
                newatoms.symbols = occs
                pop.append(newatoms)
            
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

def variable_tournament_selection(population, scores, k=6, inputseed=None):
    if inputseed != None:
        np.random.seed(inputseed)
    selection_ix = np.random.randint(len(population))
    for ix in np.random.randint(0, len(population), k-1):
        # check if better AND dissimilar (e.g. perform a tournament)
        if scores[ix] < scores[selection_ix] and scores[ix] != scores[selection_ix]:
            selection_ix = ix
    return population[selection_ix]


class Selector:
    def __init__(self,selection_style = 'tournament'):
        self.selection_style = selection_style
        self.set_selector()

    def set_selector(self):
        if self.selection_style == 'tournament':
            self.selection = tournament_selection
        if self.selection_style == 'novelty_tournament':
            self.selection = variable_tournament_selection

#MUTATION AND CROSSOVER FUNCTIONS
def crossover(p1, p2, inputseed = None, types=None, endpoint_compositions=False):
    if endpoint_compositions or types==None:
        assert len(p1) == len(p2), "parents must have the same length"
        psize = len(p1)
        if inputseed != None:
            np.random.seed(inputseed)
        cross_point = np.random.randint(1, psize-1)
        c1 = p1[:cross_point] + p2[cross_point:]
        c2 = p2[:cross_point] + p1[cross_point:]
    else:
        assert len(p1) == len(p2), "parents must have the same length"
        itr = 0
        comps_dct1 = get_comp(p1,types)
        comp_vals1 = list(comps_dct1.values())
        comps_dct2 = get_comp(p2,types)
        comp_vals2 = list(comps_dct2.values())
        while itr == 0 or any([icomp == 0.0 for icomp in comp_vals1]) or any([icomp == 0.0 for icomp in comp_vals2]):
            psize = len(p1)
            if inputseed != None:
                np.random.seed(inputseed)
            cross_point = np.random.randint(1, psize-1)
            c1 = p1[:cross_point] + p2[cross_point:]
            c2 = p2[:cross_point] + p1[cross_point:]
            comps_dct1 = get_comp(c1,types)
            comp_vals1 = list(comps_dct1.values())
            comps_dct2 = get_comp(c2,types)
            comp_vals2 = list(comps_dct2.values())
            itr += 1
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


def flip_one_atom(atoms,types,endpoint_compositions=False):
    if endpoint_compositions:
        new_atoms = atoms.copy()
        flip_ind = np.random.randint(0,len(atoms))
        flip_current = new_atoms[flip_ind].symbol
        excluded = [typ for typ in types if typ != flip_current]
        flip_to_ind = np.random.randint(0,len(excluded))
        flip_to_type = excluded[flip_to_ind]
        new_atoms[flip_ind].symbol = flip_to_type
    else:
        itr = 0
        comps_dct = get_comp(atoms,types)
        comp_vals = list(comps_dct.values())
        while itr == 0 or any([icomp == 0.0 for icomp in comp_vals]):
            new_atoms = atoms.copy()
            flip_ind = np.random.randint(0,len(atoms))
            flip_current = new_atoms[flip_ind].symbol
            excluded = [typ for typ in types if typ != flip_current]
            flip_to_ind = np.random.randint(0,len(excluded))
            flip_to_type = excluded[flip_to_ind]
            new_atoms[flip_ind].symbol = flip_to_type
            comps_dct = get_comp(new_atoms,types)
            comp_vals = list(comps_dct.values())
            itr += 1
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

def get_comp(atoms,symbols):
    comps = {symbol: 0.0 for symbol in symbols}
    counts = {symbol: 0 for symbol in symbols}
    atsymbols = [atom.symbol for atom in atoms]
    for atsymbol in atsymbols:
        counts[atsymbol] +=1
    for symbol in symbols:
        comps[symbol] = counts[symbol]/len(atoms)
    return comps

def flip_N_atoms(atoms,types,fraction=None,endpoint_compositions=False):
    if endpoint_compositions:
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
    else:
        comps_dct = get_comp(atoms,types)
        comp_vals = list(comps_dct.values())
        itr = 0
        while itr == 0 or any([icomp == 0.0 for icomp in comp_vals]):
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
            comps_dct = get_comp(new_atoms,types)
            comp_vals = list(comps_dct.values())
            itr += 1

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
