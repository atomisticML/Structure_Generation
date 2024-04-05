#!/usr/bin/env python
import sys
sys.path
import jax.numpy as np
import numpy as vnp
from jax import grad, jit, vmap
from jax import random
from functools import partial
import scipy as sp
from opt_tools import *
# get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib
import matplotlib.pyplot as plt
import lammps
import lammps.mliap
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
n_bispec=11
n_totconfig=50
#mask=[0,1,2,3,4]
mask = list(range(n_bispec))
#mask=[random.randint(0,(n_bispec-1)),random.randint(0,(n_bispec-1)),random.randint(0,(n_bispec-1)),random.randint(0,(n_bispec-1)),random.randint(0,(n_bispec-1))]
target=vnp.zeros((len(mask)))
rho0=(1/23)  #4/(3.524**3) #set a reference density
#n_atoms=random.randint(100,400)
n_atoms=5
target_comps = {'Cr':1.0}
box_rlx = False
run_temp = False
fire_min = False
line_min = True
if not fire_min and box_rlx:
    min_typ_global = 'box'
else:
    min_typ_global = 'min'

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
        self.self_entropy_grad=grad(self.self_entropy)
        self.cross_entropy_grad=grad(self.cross_entropy)
        self.V_grad=grad(self.V)
        self.K_self=self_weight
        self.K_cross=cross_weight
        self.whitening=np.identity(self.n_descriptors_keep)
        self.mode="update"
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
        dt=self.transform(d)
        self.q_count+=d.shape[0]
        #print("UPDATE: ",self.q_count)
        #update sums, covariances, and whitening matrices
        self.sum+=vnp.sum(dt,axis=0)
        self.sum_of_products+=dt.T@dt
        self.cov=self.sum_of_products/self.q_count - vnp.outer(self.sum/self.q_count, self.sum/self.q_count)
        #print("UPDATE: ",cov,self.sum/self.q_count)
        try:
            self.whitening = sp.linalg.cholesky(np.linalg.inv(self.cov))
        except ValueError:
            print('did not whiten')
            self.whitening = np.identity(self.n_descriptors_keep)
        #print("UPDATE: ",cov,self.whitening)=
        #print(list(dt))
        self.data+=list(dt)
        #print("MEAN: ", self.sum/self.q_count )
        #print("WHITENING: ", self.whitening)
        self.last_update=list(dt)
        #mean and whitening can be sent to other processes

    def self_entropy(self, descriptors):
        n_atoms=descriptors.shape[0]
        n_descriptors=descriptors.shape[1]
        self_S=0.
        d=self.transform(descriptors)
        #whitened descriptors
        wd=(self.whitening@d.T).T
        #print(whitening.shape,descriptors.shape,wd.shape)
        #Silverman's rule to set the kernel bandwidth
        #assumes that the data has been whitened
        H=vnp.identity(n_descriptors)
        H*=((n_atoms**(-1./(n_descriptors+4.)))*(4./(n_descriptors+2.)**(1./(n_descriptors+4.)) ))**2
        Hi=vnp.linalg.inv(H)
        dH=vnp.linalg.det(H)
        #covariance-multiplied descriptors
        wdh=(Hi@wd.T).T
        #this does not matter here, since it contributes only an additive constant to the entropy
        #prefactor=((2*np.pi)**(-n_descriptors/2.))*(dH**(-1./2.))
        for i in range(n_atoms):
            rho=self.density(i,wd,wdh)
            """
            rho=0.
            for j in range(n_atoms):
                if not i==j:
                    rho+=np.exp(-0.5*(wd[i]-wd[j])@(wdh[i]-wdh[j]))
            rho/=n_atoms
            """
            self_S+=-np.log(rho)

        self_S/=n_atoms
        #print(self_S)
        return self_S

    @partial(jit, static_argnums=(0,1))
    def density(self,i,wd,wdh):
        n_atoms=wd.shape[0]
        wdr=wd[i,None]-wd
        wdhr=wdh[i,None]-wdh
        rho=np.sum(np.exp(-0.5*np.sum(wdr*wdhr,axis=1)))/n_atoms

        return rho

    #this assumes that 'whitening' whitens the q data, so that q's covariance is a unit matrix
    @partial(jit, static_argnums=(0,))
    def cross_entropy(self,descriptors):
        n_atoms=descriptors.shape[0]
        n_descriptors=descriptors.shape[1]
        d=self.transform(descriptors)
        #target = vnp.array([9542.69, 44.5702, 2.49553, 39.8884, -1.30737])
        fcc_ni = vnp.array([9542.69, 44.5702, 2.49553, 39.8884, -1.30737, 1.48893, 1.77703, 13.7615, -1.17025, 0.513681, 0.764799, 0.459806, 0.87641, 34.0911, -1.39383, 1.27252, 2.42563, 0.547689, 3.454, 0.911254, 1.73656, 4.81041, 20.3444, -2.21226, 0.759079, 4.75227, 0.87024, 2.39468, 0.543481, 3.44149, 0.864251, 2.41368, 88.7406, -12.2924, 3.30387, 9.34044, 4.68315, 2.36131, 6.4768, 4.60033, 3.21413, 8.68294, 875.788, -43.9081, 25.1421, 9.73814, 13.8178, 9.23009, 18.9293, 45.946, 1624.2, 20.5024, 1.00006, 21.8708, 112.22])[:11]
        target = fcc_ni[mask]
        target=self.transform(target)
        mean=self.sum/self.q_count
        #the prefactor does not matter, as it contributes only an additive constant to the entropy
        cross_S=0.
        for i in range(n_atoms):
            x=self.whitening@(d[i]-mean)
            #x=self.whitening@(d[i])
            #PREVIOUS
            #x=self.whitening@((d[i]-target))
            #rho=np.exp(-0.5*x@x)
            #cross_S+=-np.log(rho)
            cross_S+=0.5*x@x
        cross_S/=n_atoms
        #print(cross_S)
        return cross_S

    @partial(jit, static_argnums=(0,))
    def V(self,descriptors, K_self, K_cross):
        #print("ENTROPIES: ",self.self_entropy(descriptors),self.cross_entropy(descriptors))
        #return self.K_self*self.self_entropy(descriptors)+self.K_cross*self.cross_entropy(descriptors)
        #print("IN V: ", K_self)
        #return K_self*self.self_entropy(descriptors)+K_cross*self.cross_entropy(descriptors)
        #return K_self
        return K_cross*self.cross_entropy(descriptors)+K_self*self.self_entropy(descriptors)
        #return self.self_entropy(descriptors)

    def __call__(self, elems, bispectrum, beta, energy):
        self.last_bispectrum=bispectrum.copy()
        #print("CALL ", self.mode)
        if self.mode=="run":
            #print("IN RUN")
            #print(elems.shape, bispectrum.shape, beta.shape, energy.shape)
            b=bispectrum[:,self.mask]
            #print("masked",b.shape,b)
            #print("K_self ",self.K_self, " K_cross ", self.K_cross )
            ener=self.V(b, self.K_self,self.K_cross )
            #print(ener,self.K_self,self.K_cross)
            energy[:]=0
            energy[0]=ener
            b=self.V_grad(b,self.K_self,self.K_cross)
            if not np.all(np.isfinite(b)):
                print("GRAD ERROR!")
                #print(b)

            #print(b)
            beta[:,:]=0
            beta[:,self.mask]=b
            #print(energy,b)
            #print(beta)

        if self.mode=="update":
            b=bispectrum[:,self.mask]
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
        self.lmp.commands_string(before_loading_init)
        lammps.mliap.load_model(em)
        #self.lmp.commands_string(after_loading_init)

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
box tilt large


region cell prism 0 {} 0 {} 0 {} {} {} {}
create_box 1 cell
create_atoms 1 random {} {} NULL

log log_{}.lammps
mass 1 51.9961


pair_style hybrid/overlay soft 2.0 mliap model mliappy LATER descriptor ace coupling_coefficients.yace
pair_coeff * * soft 100
pair_coeff * * mliap Cr

#variable prefactor equal ramp(0,100)
#fix 1 all adapt 1 pair soft a * * v_prefactor

thermo 10
fix nve all nve
fix lan all langevin 5000 100 1.0 48279

velocity all create 10000 4928459 dist gaussian
"""
print (generate)
def generate_cell(n_atoms, density,index):
    #shape=[5,5,5]
    shape=[10,10,10]
    xx=random.uniform(0.66*shape[0],1.33*shape[0])
    yy=random.uniform(0.66*shape[1],1.33*shape[1]+1)
    zz=random.uniform(0.66*shape[2],1.33*shape[2]+1)
    volume=xx*yy*zz
    target_volume=n_atoms/density
    scale=(target_volume/volume)**0.333333
    xx*=scale
    yy*=scale
    zz*=scale

#    dx=random.uniform(-xx,xx)*0.5
#    dy=random.uniform(-yy,yy)*0.5
#    dz=random.uniform(-zz,zz)*0.5
    dx=0.0
    dy=0.0
    dz=0.0

    seed=random.randint(1,600000)
    s=generate.format(xx,yy,zz,dx,dz,dz,n_atoms,seed,index,random.uniform(0.0,1.0))

    return s

try:
    shutil.rmtree(data_path)
except:
    pass
try:
    os.mkdir(data_path)
except:
    pass


em=EntropyModel(1,n_bispec,mask=mask)


#bootstrap the model
#g=generate_cell(n_atoms,rho0,0)
g = internal_generate_cell(0,n_atoms,template=None,desired_comps=target_comps,min_typ=min_typ_global)
sampler=EntropySampler(em,g)
sampler.update_model()

#fixes some number of atoms in the entropy minimization phase to avoid converging to perfect crystals too often.
n_fixed=0
i=0
while i <= n_totconfig:
    #rho=rho0
    print(i,"/",n_totconfig,"Using indicies :",mask)
    rho=random.uniform(0.5*rho0,1.5*rho0)
    #n_atoms=random.randint(5,50)
    #g=generate_cell(n_atoms,rho,i)
    g = internal_generate_cell(i,n_atoms,template=None,desired_comps=target_comps,min_typ=min_typ_global)
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


ad=np.array(em.data)
ad.shape
mean=(em.sum/em.q_count)
ads= ad - mean
c1=np.dot(ads.T,ads)/ads.shape[0]

#df = pd.DataFrame(ads, columns=mask)
#print(np.mean(ads,axis=0))

# print(df)

#a=pd.plotting.scatter_matrix(df)

#print(np.min(ad,axis=0))
#print(np.max(ad,axis=0))
#print(np.max(ad,axis=0)-np.min(ad,axis=0))

#c1=np.cov(ads)
#print("covariance: ",c1)

wd=(em.whitening@ads.T).T
#df = pd.DataFrame(wd, columns=mask)
#a=pd.plotting.scatter_matrix(df)

c2=np.dot(wd.T,wd)/wd.shape[0]
#print("w covariance: ",c2)
#print("w covariance: ",np.cov(wd.T))


#print(np.mean(wd,axis=0))


pickle.dump( em.data, open( "Pickled_EM-%f_%i_%f_%f.p" % (rho,n_atoms,em.K_self,em.K_cross)  , "wb" ) )
