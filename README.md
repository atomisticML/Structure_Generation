# Structure_Generation

## Requirements
FitSNAP - v x.x 
ASE - v x.x
lammps (ace mliap branch from goff)

You must add the fitsnap folder to your PYTHONPATH in order to use
the ACE descriptor setup functionalities. Otherwise, you must bring your own
ACE descriptor file (.yace).

Add the `/path/to/Structure_Generation/lib` folder to your PYTHONPATH

`export PYTHONPATH=$PYTHONPATH:/path/to/Structure_Generation/lib`

*NOTE fitsnap must also be in your pythonpath*

## Parent folder contents

`lib`   - folder containing modules : functions and classes for structure generation & optimization routines
`tools` - collecting helper scripts here
`lammps_GSQS_arbitrary` : current recommended workflow for GSQS method - uses lammps mliappy model

*deprecated/or development GSQS examples*
`genalg_GSQS` -  folder containing an important example (reproducing traditional SQS for alloys)
`genalg_GSQS_dynamic` - folder featuring the unique capability of ACE-based GSQS (ability to move move atoms off-lattice to find candidate structures for matching some target) note that this uses the 
`monte_GSQS` - similar to other GSQS example folders, but contains a simulated annealing algorithm instead of a
genetic algorithm to optimize loss functions

`other_enumeration` - other methods for structure generation & enumeration (only included for reference and comparison)
will be removed in later versions

## Definition of GSQS

For theoretical definitions and methods, see the manuscript at: (overleaf for now)

## USAGE

follow workflow and notes in `lammps_GSQS_arbitrary` 

