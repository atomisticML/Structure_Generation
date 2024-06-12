# full local minima workflow (LAMMPS) workflow in one spot

## Inputs

```
# strength of soft potential (correcting core repulsion problem in ML ACE/SNAP fingerprints)
soft_strength = 1.0
# total number of atomic configurations to generate (candidates proposed through randomly populated hermite normal form supercells)
n_totconfig=10
# max and minimum sizes for cells hermite normal form trace (equal to # of atoms in supercell)
maxcellsize=40
mincellsize=4
# you may add a mask for specific descriptors to use - here all are used 
mask = list(range(n_descs))
#target composition in dictionary form: e.g. ( W:0.5, Be:0.5 )
target_comps = {'Cr':0.7,'Fe':0.2,'Si':0.05,'V':0.05}
# extra penalty/bonus for exact matches of ACE/SNAP fingerprints with target distribution for:
# fingerprint mean within current structure compared to target distribution mean
q1_wt=5.0
# fingerprint mean within set of structures compared to target distribution mean
q1_cross_wt=5.0
# fingerprint variance within current structure compared to target distribution variance
q2_wt=0.0
# fingerprint variance within set of structures compared to target distribution variance
q2_cross_wt=0.0
# flag to perform minimization with variable cell lengths
box_rlx = False
#LAMMPS MINIMIZATION FLAGS
# flag to run at temperature first
run_temp = False
# flag to use fire style minimization (good for phase transitions & large minimizations)
fire_min = False
# flag to use linestyle quadratic minimization
line_min = True
# flag to use randomized compositions for elements in the dictionary: target_comps = {'Cr':1.0 }
randomize_comps=False
```

## Running

run with:

`python GSQS_protocol.py`
