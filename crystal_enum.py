from ase.io import read,write
import spglib as spg
from ase.lattice.cubic import *
from ase.lattice.tetragonal import *
from ase.lattice.orthorhombic import *
from ase.lattice.monoclinic import *
from ase.lattice.triclinic import *
from ase.lattice.hexagonal import *
from ase.neighborlist import *
from ase import Atoms,Atom


def optimal_bond_to_latparam(optimal_bond_length,atoms,lattice_params,tol=0.15):
	#ideal reference would be fcc Cu while finding optimal starting lattice parameters for simple cubic Cu
	tst_atoms = atoms.copy()
	natural_cutoff = max([ optimal_bond_length, np.average(natural_cutoffs(atoms))])
	atinds = [atom.index for atom in atoms]
	at_dists = {i:[] for i in atinds}
	all_dists = []
	nl = primitive_neighbor_list('ijdD',pbc=tst_atoms.pbc,positions=tst_atoms.positions ,cell=atoms.get_cell(),cutoff=natural_cutoff)
	for i,j in zip(nl[0],nl[-1]):
		at_dists[i].append(j)
		all_dists.append(j)
	if len(all_dists) == 0:
		current_bond_length = 1.2 * optimal_bond_length
	else:
		current_bond_length = np.average(all_dists)
	while current_bond_length > optimal_bond_length + tol:
		if current_bond_length > optimal_bond_length + tol:
			cell = tst_atoms.get_cell()
			cell_vec_sizes = [np.linalg.norm(v) for v in cell]
			mx = max(cell_vec_sizes)
			mx_ind = cell_vec_sizes.index(mx)
			isize = cell_vec_sizes[mx_ind]
			iratio = optimal_bond_length/isize
			istep = (  iratio * 0.2 ) #1/5 th of the size
			assert istep < 1., "check your step size %f" % istep
			new_cell = (1 - istep) * cell
			tst_atoms.set_cell(new_cell)
		elif current_bond_length < optimal_bond_length - tol:
			cell = tst_atoms.get_cell()
			cell_vec_sizes = [np.linalg.norm(v) for v in cell]
			mx = max(cell_vec_sizes)
			mx_ind = cell_vec_sizes.index(mx)
			isize = cell_vec_sizes[mx_ind]
			iratio = optimal_bond_length/isize
			istep = (  iratio * 0.2 ) #1/5 th of the size
			assert istep < 1., "check your step size %f" % istep
			new_cell = (1 + istep) * cell
			tst_atoms.set_cell(new_cell)
		#all_dists = []
		print ('istep',istep,current_bond_length, optimal_bond_length,tst_atoms.get_cell())
		nl = primitive_neighbor_list('ijdD',pbc=tst_atoms.pbc,positions=tst_atoms.positions ,cell=tst_atoms.get_cell(),cutoff=natural_cutoff)
		if len(nl[-1]) == 0:
			current_bond_length = optimal_bond_length * 1.2
		else:
			current_bond_length = np.average(nl[-1])
	
	return tst_atoms
		
def get_primitive_cell(atoms):
	cell = atoms.get_cell()
	scpos = atoms.get_scaled_positions()
	numbers = atoms.numbers
	spgcell = (cell,scpos,numbers)
	symmetry = spg.get_symmetry(spgcell, symprec=1e-5)
	lattice, scaled_positions, numbers = spg.standardize_cell(spgcell, to_primitive=True, no_idealize=False, symprec=1e-5)
	new_prim = Atoms(numbers)
	new_prim.set_cell(lattice)
	new_prim.set_scaled_positions(scaled_positions)
	new_prim.set_pbc(True)
	return new_prim

allowed_lattice_params = {
('cubic','sc'):[('a',)],
 ('cubic','bcc'):[('a',)],
 ('cubic','fcc'):[('a',)],
 ('cubic','diamond'):[('a',)],
 ('tetragonal','st'):[('a','c/a')],
 ('tetragonal','ct'):[('a','c/a')],
 ('orthorhombic','so'):[('a','b/a','c/a')],
 ('orthorhombic','baco'):[('a','b/a','c/a')],
 ('orthorhombic','fco'):[('a','b/a','c/a')],
 ('orthorhombic','boco'):[('a','b/a','c/a')],
 ('monoclinic','sm'):[('a', 'b/a', 'c/a', 'alpha')],
 ('monoclinic','bcm'):[('a', 'b/a', 'c/a', 'alpha')],
 ('triclinic','t'):[('a', 'b/a', 'c/a', 'alpha', 'beta', 'gamma')],
 ('hexagonal','h'):[('a','c/a')],
 ('hexagonal','hcp'):[('a','c/a')],
 ('hexagonal','hgr'):[('a','c/a')]
}

def lattice_func(pltup):
	valid_tups = [ ('cubic','sc'),
 ('cubic','bcc'),
 ('cubic','fcc'),
 ('cubic','diamond'),
 ('tetragonal','st'),
 ('tetragonal','ct'),
 ('orthorhombic','so'),
 ('orthorhombic','baco'),
 ('orthorhombic','fco'),
 ('orthorhombic','boco'),
 ('monoclinic','sm'),
 ('monoclinic','bcm'),
 ('triclinic','t'),
 ('hexagonal','h'),
 ('hexagonal','hcp'),
 ('hexagonal','hgr')]
	all_tup_str = ' '.join( '("%s" , "%s")'%B for B in valid_tups)
	assert pltup in valid_tups, "(%s,%s) is not a valid structure tuple, please enter one of the following: %s" % (pltup + (all_tup_str,))
	# cubic structures
	if pltup == ('cubic','sc'):
		return SimpleCubic
	elif pltup == ('cubic','bcc'):
		return BodyCenteredCubic
	elif pltup == ('cubic','fcc'):
		return FaceCenteredCubic
	elif pltup == ('cubic','diamond'):
		return Diamond
	# tetragonal structures
	elif pltup == ('tetragonal','st'):
		return SimpleTetragonal
	elif pltup == ('tetragonal','ct'):
		return CenteredTetragonal
	# orthorhombic structures
	elif pltup == ('orthorhombic','so'):
		return SimpleOrthorhombic
	elif pltup == ('orthorhombic','baco'):
		return BaseCenteredOrthorhombic
	elif pltup == ('orthorhombic','fco'):
		return FaceCenteredOrthorhombic
	elif pltup == ('orthorhombic','boco'):
		return BodyCenteredOrthorhombic
	# monoclinic structures
	elif pltup == ('monoclinic','sm'):
		return SimpleMonoclinic
	elif pltup == ('monoclinic','bcm'):
		return BaseCenteredMonoclinic
	# triclinic structure
	elif pltup == ('triclinic','t'):
		return Triclinic
	# hexagonal structures
	elif pltup == ('hexagonal','h'):
		return Hexagonal
	elif pltup == ('hexagonal','hcp'):
		return HexagonalClosedPacked
	elif pltup == ('hexagonal','hgr'):
		return Graphite


class Crystal_Enum:
	def __init__(self,base_species):
		self.a = None
		self.crystal_structures = {}
		self.base_species = base_species
		return None

	def set_lattices_per_species(self,inpdict):
		assert len(inpdict.keys()) == len(base_species), "must input a dictionary the length of the # of chemical species"
		possible_lattices = ['bcc']
		self.crystal_structures = inpdict
		for chemical, lattices_by_chemical in inpdct.items():
			self.enumerate_lattices(chemical,lattices_by_chemical)

	#def enumerate_lattices(self,chemical,lattices):


#slab = BodyCenteredCubic(directions=[[1,0,0], [0,1,0], [1,1,1]],
#			  size=(2,2,3), symbol='Cu', pbc=(1,1,0),
#			  latticeconstant=4.0)

this_func = lattice_func(('cubic','bcc'))
this_func = lattice_func(('orthorhombic','so'))
atoms = this_func( size=(1,1,1), symbol='Cu', pbc=(1,1,1), latticeconstant={'a':4.0, 'b/a':1.2, 'c/a':1.3})
prim_atoms = get_primitive_cell(atoms)
starting_atoms = optimal_bond_to_latparam(optimal_bond_length=1.5,atoms=prim_atoms,lattice_params=None,tol=0.05)
#prim_atoms = get_primitive_cell(atoms)
#print (prim_atoms)
#print  (prim_atoms.numbers)
write('atoms.vasp',starting_atoms)
