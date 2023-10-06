from yamlpace_tools.potential import *

# coupling coefficients for different bases
# number of species in the system
nspecies = 1
species = list(range(nspecies))
reference_ens = [0.]* nspecies
bonds = [b for b in itertools.product(species,repeat=2) ]
print (bonds)
rcutfaci = [4.5]*len(bonds)
lmbdai=[2.3]*len(bonds)
rcinneri = [0.0]*len(bonds)
drcinneri = [0.01]*len(bonds)

exhaust_bonds = bonds.copy()

rcutfac = []
lmbda = []
rcinner = []
drcinner = []

for bond in exhaust_bonds:
	srt_bnd = tuple(sorted(bond))
	idx = bonds.index(srt_bnd)
	rcutfac.append(rcutfaci[idx])
	lmbda.append(lmbdai[idx])
	rcinner.append(rcinneri[idx])
	drcinner.append(drcinneri[idx])
print('rcutfac =',' '.join(str(k) for k in rcutfac))
print('lambda =',' '.join(str(k) for k in lmbda))
print('rcinner =',' '.join(str(k) for k in rcinner))
print('drcinner =',' '.join(str(k) for k in drcinner))

elmap = {0:'A',1:'B',2:'C',3:'D',4:'E'}

elements = [elmap[i] for i in range(nspecies)]

# same radial and angular resolution for all #s of chemicals
ranks = [1, 2, 3, 4]
lmax =  [1, 2, 2, 1]
nmax = [4, 1, 1, 1]
lmin = 1
nradbase=4


Apot = AcePot(elements,reference_ens,ranks,nmax,lmax,nradbase,rcutfac,lmbda,rcinner,drcinner,lmin,RPI_heuristic='root_SO3_span')
# read the potential file to get expansion coefficients
Apot.write_pot('coupling_coefficients')
print ('total descs',len(Apot.nus))
