#Refactor code (GSQS_protocol.py) to get rid of redundant/cumbersome scripts
# move to a more object-oriented setup like in FS
# use a class instead of inline definitions for things done in the GRS method
# the first step is to do this for the starting structures
# to get started, look at the code in 'fitsnap3lib/solvers' (solver.py and ridge.py)
# solver.py contains the parent class for solving linear systems of equations in FS
# ridge.py is a subclass of the solver class in solver.py that does ridge regression specifically

# here, we want  a 'structure' parent class that we can use for any 'starting structure type' that we
#  may use in the GRS method
#  we will start this by making the parent structure class and then we will make 
#  a structure subclass for 'internal_generate_cell'

class Structure:
    def __init__(self)

#Make subclasses of this (like the subclasses in fitsnap) for different candidate types


    #should contain information about the type of candidate structures (crystalline, alloy, random/amorphous, etc)
