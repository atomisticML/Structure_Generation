import numpy as np
import itertools 

def get_hnfs(hnf_trs,max_tr=None, min_tr = None):

    def get_acf():
        if min_tr ==None and max_tr == None:
            acfs = []
            tracking = []
            for num in hnf_trs:
                facs = [ i for i in list(range(num+1))[1:] if num %i ==0]
                cmbs = itertools.combinations_with_replacement(facs,3)
                acf = [c for c in cmbs if np.prod(c) == num]
                for c in acf:
                    for p in itertools.permutations(c):
                        t = ','.join(str(k) for k in p)
                        if t not in tracking:# and p[1] >= num/p[0]:
                            assert p[2] == num/(p[0]*p[1]), " 'f' criterion not met"
                            acfs.append(list(p))
                            tracking.append(t)
        #else:
            #enumerate all

        return acfs,tracking
    acfs,string = get_acf()

    def get_bde(acfs,string):
        mats = []
        acf_bde ={ st:[] for st in string}
        for acf in acfs:
            bdes = []
            for b in range(acf[1]):
                for d in range(acf[2]):
                    for e in range( acf[2]):
                        bdes.append([b,d,e])
                        mat = np.array([[acf[0],0,0],
                                        [b,acf[1],0],
                                        [d,e,acf[2]]])
                        mats.append(mat)
            acf_bde[','.join(str(k) for k in acf)].append(bdes)
        return mats
    hnf_mats = get_bde(acfs,string)
    #TODO reduce by symmetry
    return hnf_mats


def limit_mats_len(hnfs,N,tol = None ):
    if tol == None:
        tol = N-1
    keeps = []
    for mat in hnfs:
        norms = np.array([np.linalg.norm(k) for k in mat])
        norms /= (N/3)
        flag = all([norm <= tol for norm in norms])
        if flag:
            keeps.append(mat)
    return keeps
N=40
full_hnf_mats = get_hnfs(hnf_trs = [N])
print ('full',len(full_hnf_mats))
hnf_mats_sub = limit_mats_len(full_hnf_mats,N,tol =0.5 )
print ('hnf_mats sub',len(hnf_mats_sub))
hnf_mats = hnf_mats_sub.copy()
for mat in hnf_mats:
    norms = np.array([np.linalg.norm(k) for k in mat])
    norms /= (N/3)
    print (norms)
print ('hnf_mats sub',len(hnf_mats_sub))
