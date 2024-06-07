from opt_tools import *
from ase.io import read,write
from ase.atoms import Atoms
#process
def set_target_struct(inp):
    if type(inp) == 'str':
        try:
            inp = read(inp)
        except:
            raise TypeError("unrecognized file type %s" % inp)
    elif type(inp) == Atoms:
    #    except 
        av,var = build_target(inp)
        
from ase.build import bulk
atoms = bulk('Cr')
#set_target_struct(atoms)

from crystal_enum import *

cell_type = get_cell_type(atoms)
print(cell_type)
