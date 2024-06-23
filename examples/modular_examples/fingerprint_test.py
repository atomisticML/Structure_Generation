
from GRS.fingerprints.ace import *

fingerprint_settings = {
'elements':['Cr','Fe'],
'ranks':[1,2,3],
'nmax':[6,2,2],
'lmax':[0,2,2],
#optional args
#'lmin':[0,0,0],
#'L_R':0,
#'M_R':0
}


fingerprint_space = ACE_ChemPos('ace',fingerprint_settings)

fingerprint_space.enumerate_labels()

print('printing fingerprint labels')
print(fingerprint_space.labels)

# write default lammps format labels
fingerprint_file_prefix = 'coupling_coefficients'
fingerprint_space.write_fingerprint_params(fingerprint_file_prefix)
