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
        norm = np.linalg.norm(mat)
        norms = np.linalg.norm(mat/norm,axis=1)
        #norms = np.array([np.linalg.norm(k) for k in mat])
        #norms /= (N)
        flag = all([norm >= tol for norm in norms])
        if flag:
            keeps.append(mat)
    return keeps

"""
def limit_mats_len(hnfs,N,tol =0.25 ):
    keeps = []
    for mat in hnfs:
        norms = np.array([np.linalg.norm(k) for k in mat])
        norms /= (N/3)
        #print (norms)
        flag = all([norm >= tol for norm in norms])
        if flag:
            keeps.append(mat)
    return keeps
"""
"""        
N=40
full_hnf_mats = get_hnfs(hnf_trs = [N])
hnf_mats_sub = limit_mats_len(full_hnf_mats,N,tol =0.16 )
hnf_mats = hnf_mats_sub.copy()
for mat in hnf_mats:
    print (mat)

#lines below can be removed if you do not want to visualize the cell
import ase
from ase.io import read,write
from ase import Atom,Atoms
from ase.build import bulk
atoms = bulk('Cu','fcc', a = 3.5)
#atoms =read('prim.cif')
from ase.visualize import view
count =1
for hnf in hnf_mats:
    det = np.linalg.det(hnf)
    cell_c = atoms.get_cell().copy()
    cell_col = cell_c.T
    hnf_cell = np.matmul(cell_col,hnf).T
    cell_pars = ase.geometry.cell_to_cellpar(hnf_cell, radians=False)
    rep = atoms *(int(det),1,1)

    new_atoms = []
    sc_mat = np.linalg.inv(hnf)
    for i,p in enumerate(atoms.get_scaled_positions()):
        print (p)
        sc_p = np.matmul(p,sc_mat)
        sc_p_add = np.matmul(p +  np.array([0.5,0.5,0.5]), sc_mat)
        new_atoms.append(Atom(atoms[i].symbol, np.matmul(sc_p , hnf_cell)))
        new_atoms.append(Atom(atoms[i].symbol, np.matmul(sc_p_add , hnf_cell)))
        
    new_atoms = Atoms(new_atoms)
    new_atoms.set_cell(hnf_cell)
    atoms.set_pbc=(True)
    renderer = write('cell_%d.pov' % count, new_atoms)
                 #**generic_projection_settings,
                 #povray_settings=povray_settings)
    renderer.render()
    #write('cell_%d.pov' % count,new_atoms,run_povray=True,**{'show_unit_cell':2})
    count +=1
"""
tmp = None
