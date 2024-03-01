import sys
from opt_tools import *
sys.path
import jax.numpy as np
import numpy as vnp
from jax import grad, jit, vmap
from jax import random
from functools import partial
import scipy as sp
# get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib
import matplotlib.pyplot as plt
import lammps
import lammps.mliap
from lammps.mliap.loader import *
import pandas as pd
import copy
import pickle
import os
import shutil
import random
#plt.rcParams['figure.figsize']=[12,8]


data_path="./StructureDump"
cross_weight=-0.001 #Entropy wrt mean of configuration. min/max just flip this sign (pos being minimize).
self_weight=100. #Entropy within each atomic configuration. min/max just flip this sign (pos being minimize). 
n_descs=11
n_totconfig=1
#mask=[0,1,2,3,4]
#mask=[random.randint(0,(n_descs-1)),random.randint(0,(n_descs-1)),random.randint(0,(n_descs-1)),random.randint(0,(n_descs-1)),random.randint(0,(n_descs-1))]
mask = list(range(n_descs))
target=vnp.zeros((len(mask)))
rho0=4/(3.524**3) #set a reference density
#n_atoms=random.randint(100,400)
n_atoms=32

descriptors_target_1 = np.load('target_descriptors.npy')

class EntropyModel:
    def __init__(self, n_elements, n_descriptors_tot, mask):
        self.mask=mask
        self.n_descriptors=n_descriptors_tot
        self.n_descriptors_keep=len(mask)
        self.n_elements=n_elements
        self.n_params=1
        self.sum_of_products=vnp.zeros((self.n_descriptors_keep,self.n_descriptors_keep))
        self.sum=vnp.zeros((self.n_descriptors_keep,))
        self.q_count=0
        self.first_moment_grad=grad(self.first_moment)
        #self.cross_entropy_grad=grad(self.cross_entropy)
        self.V_grad=grad(self.V)
        self.K_self=self_weight
        self.K_cross=cross_weight
        self.whitening=np.identity(self.n_descriptors_keep)
        self.mode=  "update"
        self.data=[]

    def set_mode_update(self):
        self.mode="update"

    def set_mode_run(self):
        self.mode="run"

    @partial(jit, static_argnums=(0,))
    def transform(self,d):
        #return np.log(d)
        return d

    def update(self,d):
        #print("UPDATE: ")
        #print(d)
        #transform the raw descriptors
        #dt=self.transform(d)
        #self.q_count+=d.shape[0]
        #print("UPDATE: ",self.q_count)
        #update sums, covariances, and whitening matrices
        #self.sum+=vnp.sum(dt,axis=0)
        #self.sum_of_products+=dt.T@dt
        #self.cov=self.sum_of_products/self.q_count - vnp.outer(self.sum/self.q_count, self.sum/self.q_count)
        #print("UPDATE: ",cov,self.sum/self.q_count)
        #self.whitening = sp.linalg.cholesky(np.linalg.inv(self.cov))
        #print("UPDATE: ",cov,self.whitening)=
        #print(list(dt))
        #self.data+=list(dt)
        #print("MEAN: ", self.sum/self.q_count )
        #print("WHITENING: ", self.whitening)
        #self.last_update=list(dt)
        #mean and whitening can be sent to other processes
        self.set_mode_run()

    def first_moment(self, descriptors):
        n_atoms=descriptors.shape[0]
        n_descriptors=descriptors.shape[1]
        avgs = np.average(descriptors,axis=0)
        abs_diffs_1 = [np.abs(ii - kk) for ii,kk in zip(avgs, descriptors_target_1)]
        abs_diffs_1 = np.array(abs_diffs_1)
        tst_residual_1 = np.sum(abs_diffs_1)
        print (tst_residual_1)
        return tst_residual_1

    @partial(jit, static_argnums=(0,1))
    def density(self,i,wd,wdh):
        n_atoms=wd.shape[0]
        wdr=wd[i,None]-wd
        wdhr=wdh[i,None]-wdh
        rho=np.sum(np.exp(-0.5*np.sum(wdr*wdhr,axis=1)))/n_atoms

        return rho

    @partial(jit, static_argnums=(0))
    def V(self,descriptors):
        return self.first_moment(descriptors)

    def __call__(self, elems, descriptors, beta, energy):
        self.last_descriptors=descriptors.copy()
        #print("CALL ", self.mode)
        if self.mode=="run":
            #print("IN RUN")
            #print(elems.shape, descriptors.shape, beta.shape, energy.shape)
            b=descriptors[:,self.mask]
            #print("masked",b.shape,b)
            #print("K_self ",self.K_self, " K_cross ", self.K_cross )
            ener=self.V(b)
            #print(ener,self.K_self,self.K_cross)
            energy[:]=0
            energy[0]=ener
            b=self.V_grad(b)
            if not np.all(np.isfinite(b)):
                print("GRAD ERROR!")
                #print(b)

            #print(b)
            beta[:,:]=0
            beta[:,self.mask]=b
            #print(energy,b)
            #print(beta)

        if self.mode=="update":
            b=descriptors[:,self.mask]
            self.update(b)

class EntropySampler:
    def __init__(self, model, before_loading_init):
        self.model=model
        #create an instance of LAMMPS
        self.lmp = lammps.lammps(cmdargs=['-screen','none'])
        # Before defining the pair style, one must do the following:
        lammps.mliap.activate_mliappy(self.lmp)
        # Otherwise, when running lammps in library mode,
        # you will get an error:
        # "ERROR: Loading MLIAPPY coupling module failure."
        # Setup the simulation and declare an empty model
        # by specifying model filename as "LATER"
        print(before_loading_init)
        self.lmp.commands_string(before_loading_init)
        print (em,model)
        lammps.mliap.load_model(em)
        #self.lmp.commands_string(after_loading_init)
        """
        self.model=model
        #create an instance of LAMMPS
        self.lmp = lammps.lammps(cmdargs=['-screen','none'])
        # Before defining the pair style, one must do the following:
        lammps.mliap.loader.activate_mliappy(self.lmp)
        print ('before mliappy')
        print (self.model,model)
        activate_mliappy(self.lmp)
        print ('after mliappy')
        # Otherwise, when running lammps in library mode,
        # you will get an error:
        # "ERROR: Loading MLIAPPY coupling module failure."
        # Setup the simulation and declare an empty model
        # by specifying model filename as "LATER"
        print (before_loading_init)
        self.lmp.commands_string(before_loading_init)
        print ('before loading model')
        #lammps.mliap.load_model(self.model)
        lammps.mliap.load_model(em)
        #lammps.mliap.load_model(model)
        #self.lmp.commands_string(after_loading_init)
        """
    def update_model(self):
        self.model.set_mode_update()
        self.lmp.commands_string("variable etot equal etotal")
        self.lmp.commands_string("variable ptot equal press")
        self.lmp.commands_string("variable pairp equal epair")
        self.lmp.commands_string("variable numat equal atoms")
        self.lmp.commands_string("run 0")
        self.lmp.commands_string("print \"${etot} ${pairp} ${ptot} ${numat} \" append Summary.dat screen no")
        self.model.set_mode_run()

    def run(self,cmd=None):
        if cmd==None:
            self.lmp.commands_string("run 0")
        else:
            self.lmp.commands_string(cmd)

generate=\
"""
units           metal
boundary        p p p

#lattice         fcc 3.524
#region          myreg block 0 2 0 2 0 2
#create_box      1 myreg
#create_atoms    1 box
#displace_atoms  all random 0.1 0.1 0.1 123456


read_data  {}

log log_{}.lammps
mass 1 58.6934


pair_style hybrid/scaled {} soft 2.0 1.0 LATER mliap model mliappy LATER descriptor pace coupling_coefficients.yace 1 0

pair_coeff * * soft 100
pair_coeff * * mliap Ni


#variable prefactor equal ramp(0,100)
#fix 1 all adapt 1 pair soft a * * v_prefactor

thermo 10
fix nve all nve
fix lan all langevin 5000 100 1.0 48279

velocity all create 10000 4928459 dist gaussian
"""
try:
    shutil.rmtree(data_path)
except:
    pass
try:
    os.mkdir(data_path)
except:
    pass
lmpi = lammps.lammps(cmdargs=['-screen','none'])
activate_mliappy(lmpi)
em=EntropyModel(1,n_descs,mask=mask)
#bootstrap the model
#g=gen_random(0)
g=internal_generate_cell(0)
"""
from ase.build import bulk
template = bulk('Ni',cubic=True)*(1,2,2)
all_species = template.symbols
cellg = template.get_cell()
rnd = starting_generation(1,all_species,cellg,typ='ase')
g = rnd[0]
"""
print ('model',em)
sampler=EntropySampler(em,g)
print ('after sampler')
sampler.update_model()
#print ('after update')

#fixes some number of atoms in the entropy minimization phase to avoid converging to perfect crystals too often.
n_fixed=0
i=1
while i <= n_totconfig:
    #rho=rho0
    print(i,"/",n_totconfig,"Using indicies :",mask)
    rho=random.uniform(0.5*rho0,1.5*rho0)
    n_atoms=random.randint(5,50)
    g = internal_generate_cell(i)
    #g=gen_random(i)
    #g = starting_generation(1,all_species,cellg,typ='ase')[0]
    sampler=EntropySampler(em,g)
    #em.K_cross=-i/2.
    em.K_cross=cross_weight #Entropy within the database.
    em.K_self=self_weight #to minimize/maximize just flip this sign (negative being maximize). Entropy within each atomic configuration.
    #em.K_self=-0.
    #sampler.run("minimize 1e-2 1e-2 1000 1000")
    #sampler.run("run 10000")
#     sampler.run("compute b all sna/atom 4.1 0.99363 8 1.0 0.50 bzeroflag 0")
#     sampler.run("write_dump all custom %s/sample.bi_spectrum.%i.dump id type x y z c_bisp[*]" % (data_path,i) )
    sampler.run("minimize 1e-8 1e-8 1000 1000")
    sampler.run("write_data %s/sample.%i.dat " % (data_path,i) )
    #em.K_cross*=-1. #Entropy within the database.
    #em.K_self*=-1. #to minimize/maximize just flip this sign (negative being maximize). Entropy within each atomic configuration.
    #sampler.run("minimize 1e-8 1e-8 1000 1000")
    #sampler.run("write_data %s/sample.max_self.%i.dat " % (data_path,i) )
    sampler.update_model()

    i+=1


#ad=np.array(em.data)
#ad.shape
#mean=(em.sum/em.q_count)
#ads= ad - mean
#c1=np.dot(ads.T,ads)/ads.shape[0]

#df = pd.DataFrame(ads, columns=mask)
#print(np.mean(ads,axis=0))

# print(df)

#a=pd.plotting.scatter_matrix(df)

#print(np.min(ad,axis=0))
#print(np.max(ad,axis=0))
#print(np.max(ad,axis=0)-np.min(ad,axis=0))

#c1=np.cov(ads)
#print("covariance: ",c1)

#wd=(em.whitening@ads.T).T
#df = pd.DataFrame(wd, columns=mask)
#a=pd.plotting.scatter_matrix(df)

#c2=np.dot(wd.T,wd)/wd.shape[0]
#print("w covariance: ",c2)
#print("w covariance: ",np.cov(wd.T))


#print(np.mean(wd,axis=0))


#pickle.dump( em.data, open( "Pickled_EM-%f_%i_%f_%f.p" % (rho,n_atoms,em.K_self,em.K_cross)  , "wb" ) )
