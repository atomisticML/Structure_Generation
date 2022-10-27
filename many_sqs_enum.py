from ase import Atom
import itertools
from ase.build import bulk
from ase.io import read,write
from icet import ClusterSpace
import numpy as np
from icet.tools.structure_generation import (generate_sqs,
                                             generate_sqs_from_supercells,
                                             generate_sqs_by_enumeration,
                                             generate_target_structure)

from icet.input_output.logging_tools import set_log_config
set_log_config(level='INFO')

base_species = ["Mo", "Nb", "Ta", "Ti", "W"]
primitive_structure = bulk('W',cubic=False)
cs = ClusterSpace(primitive_structure, [6., 4.0], base_species)
print (cs)
def gen_sqs(ind,target_concentrations,cs=cs):
	primitive_structure = bulk('W',cubic=False)
	cubic_structure = bulk('W',cubic=True)
	this_structure=cubic_structure*(5,5,5)
	base_species = ["Mo", "Nb", "Ta", "Ti", "W"]
	#HEAVY tungsten
	max_len = len(this_structure)
	#conc_rest = {elm: 1/len(base_species) for elm in base_species}
	#target_concentrations = conc_rest.copy()
	print(target_concentrations)
	#sqs = generate_sqs(cluster_space=cs,
	#try:
	sqs = generate_sqs_from_supercells(cluster_space=cs,
			   #max_size=10,
			   n_steps=max_len * 50,
			   #n_steps = 3000,
			   supercells= [this_structure],# cubic_structure*(2,2,2)],
			   target_concentrations=target_concentrations)
	write('sqs_%d.vasp' % ind,sqs)
	write('sqs_%d.cif' %ind,sqs)
	print('Cluster vector of generated structure:', cs.get_cluster_vector(sqs))
	#except ValueError:
	#	pass
	#del cs

chems = ['Mo','Nb','Ta','Ti','W']
target_concentrations_many = {
#1:{'Mo': 0.2, 'Nb': 0.2, 'Ta': 0.2, 'Ti': 0.2, 'W': 0.2},
2:{'Mo': 0.2, 'Nb': 0.2, 'Ta': 0.16, 'Ti': 0.16, 'W': 0.28},
3:{'Mo': 0.2, 'Nb': 0.16, 'Ta': 0.2, 'Ti': 0.16, 'W': 0.28},
4:{'Mo': 0.16, 'Nb': 0.2, 'Ta': 0.2, 'Ti': 0.16, 'W': 0.28},
5:{'Mo': 0.16, 'Nb': 0.2, 'Ta': 0.16, 'Ti': 0.2, 'W': 0.28},
6:{'Mo': 0.16, 'Nb': 0.16, 'Ta': 0.2, 'Ti': 0.2, 'W': 0.28},
}
for tind,tconc in target_concentrations_many.items():
	gen_sqs(tind,tconc)	
