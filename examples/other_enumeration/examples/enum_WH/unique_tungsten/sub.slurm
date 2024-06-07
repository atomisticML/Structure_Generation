#!/bin/bash
#SBATCH --job-name=symscfs
#SBATCH --account=fy150075
#SBATCH --partition=short,batch
#SBATCH --time=3:59:00
#SBATCH --nodes=1
#SBATCH --output=out_subrun

### End of SLURM options ###
source ~/.bashrc
module purge
module load cde/v3/devpack
module unload cde/v3/anaconda3/2021.05
export ESPRESSO_PSEUDO=/ascldap/users/jmgoff/projects/software/automation/qe_pseudos
export PATH=$PATH:/ascldap/users/jmgoff/projects/software/qe-7.1/build/bin
export PYTHONPATH=$PYTHONPATH:/ascldap/users/jmgoff/projects/software/automation/tools

nodes=$SLURM_JOB_NUM_NODES    # Number of nodes 
cores=16 # Number MPI processes to run on each node (a.k.a. PPN)

python subrun_sym_scfs.py
