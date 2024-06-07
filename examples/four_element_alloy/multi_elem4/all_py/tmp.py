from hnf import *
from ase.build import bulk
from ase import Atoms,Atom
def limit_mats_len(hnfs,N,tol = None ):
    if tol == None:
        tol = N-1
    keeps = []
    for mat in hnfs:
        norms = np.array([np.linalg.norm(k) for k in mat])
        norms /= (N/3)
        print (norms)
        flag = all([norm <= tol for norm in norms])
        if flag:
            keeps.append(mat)
    return keeps

desired_comps = {'Cr':1.0}
desired_size=4
chems = list(desired_comps.keys())
template = Atoms([chems[0]]*desired_size)
atoms_base = bulk(chems[0])
vol_base = np.dot(np.cross(atoms_base.get_cell()[0],atoms_base.get_cell()[1]),atoms_base.get_cell()[2])
a_simp = vol_base**(1/3)
cells_all = get_hnfs(hnf_trs=[desired_size])
#hnf_mats_sub = limit_mats_len(full_hnf_mats,N,tol =0.16 )
cells= limit_mats_len(cells_all,desired_size,tol = 2 )
cell = a_simp * cells[-2]
print(cells)
print(cell, np.dot(np.cross(cell[0],cell[1]),cell[2] )/desired_size, a_simp **3 )
