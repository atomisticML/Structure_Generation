from yamlpace_tools.potential import *

elements = ['Ca','Mg']
reference_ens = [0.,0.]
ranks = [1,2,3,4]
lmax = [1,2,2,1]
ldict = {1:1,2:2,3:2,4:1}
lrng = [1]
L_R = 0
M_R= 0

try:
    with open('wig_LR_%d_r4_lmax%d.pickle' %(L_R,max(lmax)),'rb') as handle:
        ccs = pickle.load(handle)
except FileNotFoundError:
    ccs = get_wig_coupling(ldict,L_R=L_R)
    store_generalized(ccs, coupling_type='wig',L_R=L_R)



nu_per_rank = {1:[], 2:[],3:[],4:[]}

rcutfac = [5.5, 5.75, 5.75, 6]
lmbda=[k*0.05 for k in rcutfac]


rcinner = [0.0] * len(lmbda)
drcinner = [0.0]* len(lmbda)

elements = ['Ca','Mg']
nmax = [2,2,2,1]
lmin = [0,0,0,1]
nradbase=max(nmax)
#permutation symmetry adapted ACE labels
Apot = AcePot(elements,reference_ens,ranks,nmax,lmax,nradbase=nradbase,rcut=rcutfac,lmbda=lmbda,rcutinner=rcinner,drcutinner=drcinner,lmin = lmin,**{'ccs':ccs[M_R]})


Apot.write_pot('coupling_coefficients')
