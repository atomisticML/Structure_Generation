# GSQS workflow with ACE ML-IAP

## Requirements
FitSNAP - v x.x 
ASE - v x.x
lammps (ace mliap branch from goff)

You must add the fitsnap folder to your PYTHONPATH in order to use
the ACE descriptor setup functionalities. Otherwise, you must bring your own
ACE descriptor file (.yace).

Add the `/path/to/StructureGeneration/lib` folder to your PYTHONPATH

`export PYTHONPATH=$PYTHONPATH:/path/to/StructureGeneration/lib`


## Workflow and Usage

### Step 1: system and ACE descriptor initialization

First copy the `setup_lammps_descs.py` file from the `/path/to/StructureGeneration/tools`
to your current directory.

Edit the elements, nmax, lmax, and other hyperparameters to reasonable values. Choosing
elements to a list of all elements in your target systems. There are some tools to get
ballparks for these numbers in FitSNAP/tools. After this is adjusted, run the script.

### Step 2: Set target distribution
The most straightforward usage is to provide an atomistic structure file 
(any ASE-compatible format) to use as the target. Simply run the get_target.py script
and provide the name of the target file as a command line argument. In this example
the target structure is generated with ASE by running the build_target.py script. This
produces a `supercell_target.cif` for which we must obtain the ACE descriptor
distribution for. The target ACE descriptor distribution is obtained and saved by:

`python get_target.py supercell_target.cif`

Alternatively, one can define the target distribution with a saved numpy files
for both the average descriptor values and the variances manually:

`target_descriptors.npy`

`target_var_descriptors.npy`


### Set parameters in GSQS_protocol.py and run

If relevant, adjust the minimization and other parameters in the `GSQS_protocol.py` script.
In that script, you may choose the # of structures to generate, min/max structure sizes,
and more.

provided the previous steps produced your target distributions correctly, you may now run 
the LAMMPS mliappy driver for GSQS. Run with:

`python GSQS_protocol.py`

Structures will be stored in `StructureDump` with a summary of the mliappy potential per
candidate in Summary.dat. The lammps data files in `StructureDump` may be converted to
cif, POSCAR, and other formats using the `cif_all.py` script (defaults to cif format).
Simply copy the cif_all script into the StructureDump folder and run it to obtain all of
the cif files.
