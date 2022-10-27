from ase import Atom
from ase.build import bulk
from ase.io import read,write
from icet import ClusterSpace
from icet.tools.structure_generation import (generate_sqs,
                                             generate_sqs_from_supercells,
                                             generate_sqs_by_enumeration,
                                             generate_target_structure)

from icet.input_output.logging_tools import set_log_config
set_log_config(level='INFO')

primitive_structure = bulk('W')
cubic_structure = bulk('W',cubic=True)
base_species = ["Mo", "Nb", "Ta", "Ti", "W"]
#HEAVY tungsten
max_len = 27
target_concentrations = {'Mo': 5/27, 'Nb': 5/27, 'Ta': 5/27, 'Ti': 5/27, 'W': 7/27}
#conc_rest = {elm: 1/len(base_species) for elm in base_species}
cs = ClusterSpace(primitive_structure, [6.0, 3.0], base_species)
print (cs)
#target_concentrations = conc_rest.copy()
print(target_concentrations)
#sqs = generate_sqs(cluster_space=cs,
sqs = generate_sqs_from_supercells(cluster_space=cs,
                   #max_size=10,
                   n_steps=max_len * 800,
                   supercells= [primitive_structure*(3,3,3)],# cubic_structure*(2,2,2)],
                   target_concentrations=target_concentrations)
write('sqs_1.vasp',sqs)
write('sqs_1.cif',sqs)
print('Cluster vector of generated structure:', cs.get_cluster_vector(sqs))

target_concentrations = {'Mo': 5/27, 'Nb': 5/27, 'Ta': 5/27, 'Ti': 5/27, 'W': 7/27}
