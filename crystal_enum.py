from ase.io import read,write
import spglib as spg
from ase.lattice.cubic import BodyCenteredCubic
from ase import Atoms,Atom

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

	def enumerate_lattices(self,chemical,lattices):
		
		



slab = BodyCenteredCubic(directions=[[1,0,0], [0,1,0], [1,1,1]],
			  size=(2,2,3), symbol='Cu', pbc=(1,1,0),
			  latticeconstant=4.0)

atoms = BodyCenteredCubic( size=(1,1,1), symbol='Cu', pbc=(1,1,1), latticeconstant=4.0)
print (atoms)
print  (atoms.numbers)
prim = get_primitive_cell(atoms)
write('atoms.vasp',prim)


"""
lattice.cubic

    SimpleCubic

    FaceCenteredCubic

    BodyCenteredCubic

    Diamond (*)

lattice.tetragonal

    SimpleTetragonal

    CenteredTetragonal

lattice.orthorhombic

    SimpleOrthorhombic

    BaseCenteredOrthorhombic

    FaceCenteredOrthorhombic

    BodyCenteredOrthorhombic

lattice.monoclinic

    SimpleMonoclinic

    BaseCenteredMonoclinic

lattice.triclinic

    Triclinic

lattice.hexagonal

    Hexagonal

    HexagonalClosedPacked (*)

    Graphite (*)

"""


"""
Cubic
	

a
	

‘a’

Tetragonal
	

(a, c)
	

‘a’, ‘c’ or ‘c/a’

Orthorhombic
	

(a, b, c)
	

‘a’, ‘b’ or ‘b/a’, ‘c’ or ‘c/a’

Triclinic
	

(a, b, c,
, ,

)
	

‘a’, ‘b’ or ‘b/a’, ‘c’ or ‘c/a’, ‘alpha’, ‘beta’, ‘gamma’

Monoclinic
	

(a, b, c, alpha)
	

‘a’, ‘b’ or ‘b/a’, ‘c’ or ‘c/a’, ‘alpha’

Hexagonal
	

(a, c)
	

‘a’, ‘c’ or ‘c/a’
"""
