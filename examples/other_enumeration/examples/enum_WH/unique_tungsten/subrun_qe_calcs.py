import glob
import os
from subprocess import call

struct_files = sorted(glob.glob('*.cif'))
for fileindex,sfile in enumerate(struct_files):
    print ('running',sfile)
    prefix = sfile.split('.')[0]
    dirname = 'run_%03d' % fileindex
    os.chdir(dirname)
    if not os.path.isdir('out'):
        call('mpirun -np 16 cp.x -in %s_md.in > %s_md.out ' % (prefix,prefix),shell = True)
        #call('rm ./out/*/*wfc*', shell = True)
        #call('rm ./out/*/*upf', shell = True)
        #call('rm ./out/*/*dat', shell = True)
    else:
        pass
    os.chdir('../')

