# Structure_Generation

This branch is a modular refactor of the original GRS/GSQS codes.
See the GRS folder for the Classes, Methods, and Functions to be
included in this refactor. This branc is under active development
and versions are in flux.

## Requirements
FitSNAP - (post library update) 
ASE - v 3.22.1
lammps - LAMMPS (post mliap update)

You must add the fitsnap folder to your PYTHONPATH in order to use
the ACE descriptor setup functionalities. Otherwise, you must bring your own
ACE descriptor file (.yace).

Add the `/path/to/Structure_Generation/lib` folder to your PYTHONPATH

`export PYTHONPATH=$PYTHONPATH:/path/to/Structure_Generation/lib`

*NOTE fitsnap must also be in your pythonpath*

## Parent folder contents

`GRS` - generalized representative structure library

`lib`   - folder containing helper functions and classes for structure generation & optimization routines
`tools` - some extra helper scripts here - will likely be resorted 

`examples` folder containing examples for local (lammps) and global optimization (genetic algorithm, simulated annealing) - See further details in the README of the examples folder


## Definition of GSQS

For theoretical definitions and methods, see the manuscript at: *link_here*
*TODO* update with software-relevant definitions

## USAGE

follow procedure in `examples/simple_test` first. 
