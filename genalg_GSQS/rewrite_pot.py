from yamlpace_tools.potential import *

# potential parameters. these should match the FitSNAP input
reference_ens = [0.,0.]
# bonds:
bonds=  [(0, 0), (0, 1), (1, 0), (1, 1)]
#  [(H, H), (H, N), (H, O), (N, N), (N, O), (O, O)]

rcutfaci = [4.5]* 4 #[5.0, 5.5, 5.7, 4.4, 5.7, 5.5]
lmbdai = [2.3]*4
rcinneri = [0.0]*4
drcinneri = [0.01]*4

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

elements=["Ca","Mg"]
ranks = [1, 2, 3, 4]
lmax =  [1, 2, 2, 1]
nmax = [4, 1, 1, 1]
lmin = 1
nradbase=4


Apot = AcePot(elements,reference_ens,ranks,nmax,lmax,nradbase,rcutfac,lmbda,rcinner,drcinner,lmin,RPI_heuristic='root_SO3_span')
# read the potential file to get expansion coefficients
Apot.write_pot('coupling_coefficients')
