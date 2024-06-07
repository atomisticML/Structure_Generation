#NOTE you must have FitSNAP in your python path to use this
from fitsnap3lib.lib.sym_ACE.yamlpace_tools.potential import *
from lib.default_ACE_settings import *
#flag to use generalized CG coefficients instead of generalized Wigner symbols
# FitSNAP default is 'False'.
cg_flag = True

# define the ACE descriptor set
elements = ['Cr','Fe','Si','V']
# reference value for 0th order ace descriptor 
reference_ens = [0.]*len(elements)
# ranks of descriptors
ranks = [1,2,3]
# angular character of fingerprints
lmax = [0,2,2]
# minimum l per rank (for trimming descriptor list)
lmin = [0,0,1]
# radial character of fingerprints
nmax = [3,2,1]
nradbase=max(nmax)

# radial function hyperparameters are defined per bond type:
bonds = [p for p in itertools.product(elements,elements)]
print('bonds',bonds)
#rcutfac=[4.76, 4.76, 4.75, 4.675, 4.76, 4.76, 4.75, 4.675, 4.75, 4.75, 5.4, 4.9, 4.675, 4.675, 4.9, 4.59]
#lmbda=[0.238, 0.238, 0.238, 0.234, 0.238, 0.238, 0.238, 0.234, 0.238, 0.238, 0.27, 0.245, 0.234, 0.234, 0.245, 0.23]
rc_range,rc_default,lmb_default,rcin_default = get_default_settings(elements,nshell=2,return_range=True,apply_shift=False)
rcutfac = [float(k) for k in rc_default.split()[2:]]
lmbda = [float(k) for k in lmb_default.split()[2:]]
#you typically don't want an inner cutoff for structure generation
rcinner = [0.0] * len(bonds)
drcinner = [0.0]* len(bonds)


#-------------------------------------------------------------------------------------
ldict = {rank:lmax[i] for i,rank in enumerate(ranks)}
L_R = 0
M_R= 0

istr = '%d'
lstr = ''
rankstr = ''
for i,val in enumerate(ranks):
    lstr += istr
    rankstr += istr
lstr = lstr % tuple(lmax)
rankstr = rankstr % tuple(ranks)

if cg_flag:
    from fitsnap3lib.lib.sym_ACE.clebsch_couple import *
    #from clebsch_couple import *
    try:
        with open('cg_LR_%d_r%s_lmax%s.pickle' %(L_R,rankstr,lstr),'rb') as handle:
            ccs = pickle.load(handle)
    except FileNotFoundError:
        ccs = get_cg_coupling(ldict,L_R=L_R)
        #print (ccs)
        #store them for later so they don't need to be recalculated
        store_generalized(ccs, coupling_type='cg',L_R=L_R)
else:
    from fitsnap3lib.lib.sym_ACE.wigner_couple import *
    #from wigner_couple import *
    try:
        with open('wig_LR_%d_r%s_lmax%s.pickle' %(L_R,rankstr,lstr),'rb') as handle:
            ccs = pickle.load(handle)
    except FileNotFoundError:
        ccs = get_wig_coupling(ldict,L_R=L_R)
        store_generalized(ccs, coupling_type='wig',L_R=L_R)
"""
# radial function hyperparameters are defined per bond type:
bonds = [p for p in itertools.product(elements,elements)]
print('bonds',bonds)
# you typically want very small radial lambdas - increasing lambda biases your sampling of short distances
# increasing lambda too much for structure generation will probably lead to numerical instability
#rcutfac=[4.76, 4.76, 4.75, 4.675, 4.76, 4.76, 4.75, 4.675, 4.75, 4.75, 5.4, 4.9, 4.675, 4.675, 4.9, 4.59]
#lmbda=[0.238, 0.238, 0.238, 0.234, 0.238, 0.238, 0.238, 0.234, 0.238, 0.238, 0.27, 0.245, 0.234, 0.234, 0.245, 0.23]


rc_range,rc_default,lmb_default,rcin_default = get_default_settings(elements,nshell=2,return_range=True,apply_shift=False)
rcutfac = [float(k) for k in rc_default.split()[2:]]
lmbda = [float(k) for k in lmb_default.split()[2:]]
#you typically don't want an inner cutoff for structure generation
rcinner = [0.0] * len(bonds)
drcinner = [0.0]* len(bonds)
"""
Apot = AcePot(elements,reference_ens,ranks,nmax,lmax,nradbase,rcutfac,lmbda,rcinner,drcinner,lmin=lmin, **{'ccs':ccs[M_R]})
Apot.write_pot('coupling_coefficients')

