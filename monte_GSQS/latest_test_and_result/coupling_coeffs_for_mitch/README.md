Coupling coefficients files for lammps compute:

You must specify the element types in the coupling_coefficients.yace file before running
the lammps compute. Specify the element types in the first line -dummy elements A-E are used 
by default:
```
elements: [A, B]
changed to the following for an Au/Cu system:
elements: [Au, Cu]
These should be entered in alphabetical order.
```

Using the lammps compute:
#              compute id , group, compute name (pace), coupling coefficient file name, bikflag, dgradflag
compute        desc all pace coupling_coefficients.yace 1 0 

The descriptor values, averages, etc. may easily be obtained using the python script in example
titled 'get_array.py'.



