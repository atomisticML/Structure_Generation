from fitsnap3lib.lib.sym_ACE.yamlpace_tools.potential import *
#from yamlpace_tools.potential import *

#flag to use generalized CG coefficients instead of generalized Wigner symbols
# FitSNAP default is 'False'.
cg_flag = True

elements = ['H','O']
reference_ens = [0.,0.]
ranks = [1,2,3,4]
lmax = [0,2,2,1]
# minimum l per rank
lmin = [0,0,1,1]
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

ranks = [1,2,3,4]
# radial function hyperparameters
rcutfac = [4.0,4.8,4.8,4.2]
lmbda=[0.3,0.5,0.5,0.4]
rcinner = [0.0] * 4
drcinner = [0.0]* 4
nmax = [6,2,2,1]
nradbase=max(nmax)

#can play around with manually selected descriptor labels
"""
nus = [ [], [], [], ['0_0,0,0,0,1,1,1,1,1,1,1,1_2-2'] ]

for nu in nus_f:
    mu0,mu,n,l,inter = get_mu_n_l(nu,return_L=True)
    rank = get_mu_nu_rank(nu)
    nu_per_rank[rank].append(nu)

nus = list(nu_per_rank.values())
Apot = AcePot(elements,reference_ens,ranks,nmax,lmax,nradbase,rcutfac,lmbda,rcinner,drcinner,lmin=lmin, **{'input_nus':nus,'ccs':ccs[M_R]})
Apot.write_pot('coupling_coefficients')
"""
#otherwise, enumerate all functions given the parameters above
Apot = AcePot(elements,reference_ens,ranks,nmax,lmax,nradbase,rcutfac,lmbda,rcinner,drcinner,lmin=lmin, **{'ccs':ccs[M_R]})
Apot.write_pot('coupling_coefficients')

