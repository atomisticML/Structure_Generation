from ase import Atoms
from ase.data import atomic_numbers
from ase.ga.data import PrepareDB
from ase.ga.startgenerator import StartGenerator
from ase.ga.utilities import CellBounds, closest_distances_generator

symbols =['Cr']

# Generate a dictionary with the closest allowed interatomic distances
#Z = atomic_numbers['Ag']
#nums = [Z]
nums = [atomic_numbers[sym] for sym in symbols]
blmin = closest_distances_generator(atom_numbers=nums,
                                    ratio_of_covalent_radii=0.5)
print(blmin)
