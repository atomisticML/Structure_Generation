def build_target(start):
    from ase.io import read,write
    from ase import Atoms,Atom
    from ase.ga.utilities import closest_distances_generator, CellBounds
    from ase.ga.startgenerator import StartGenerator
    from ase.data import atomic_numbers, atomic_names, atomic_masses, covalent_radii
    #from __future__ import print_function
    import sys,os
    import ctypes
    import numpy as np
    from lammps import lammps, LMP_TYPE_ARRAY, LMP_STYLE_GLOBAL

    # get mpi settings from lammps
    def run_struct(atoms,fname,maxcut=6.0):
            
        lmp = lammps()
        me = lmp.extract_setting("world_rank")
        nprocs = lmp.extract_setting("world_size")


        cmds = ["-screen", "none", "-log", "none"]
        lmp = lammps(cmdargs = cmds)

        #def set_atoms(atoms,atid=0):
        #    write('iter_%d.data' % atid,atoms,format='lammps-data')
        #    lmp.command('read_data  iter_%d.data' % atid )
        #    lmp.command('mass  1 180.94788')
        #    lmp.command(f"run {nsteps}")

        def run_lammps(dgradflag):

            # simulation settings
            fname = file_prefix
            lmp.command("clear")
            lmp.command("info all out log")
            lmp.command('units  metal')
            lmp.command('atom_style  atomic')
            lmp.command("boundary    p p p")
            lmp.command("atom_modify    map hash")
            lmp.command('neighbor  2.3 bin')
            # boundary
            lmp.command('boundary  p p p')
            # read atoms
            lmp.command('read_data  %s.data' % fname )
            utypes = []
            for atom in atoms:
                if atom.symbol not in utypes:
                    utypes.append(atom.symbol)
            for ind,typ in enumerate(utypes):
                number = atomic_numbers[typ]
                mass = atomic_masses[number]
                lmp.command('mass   %d %f' % (ind+1,mass))

            lmp.command("pair_style     zero %f" % maxcut)
            lmp.command(f"pair_coeff     * *")

            if dgradflag:
                lmp.command(f"compute     pace all pace coupling_coefficients.yace 1 1")
            else:
                lmp.command(f"compute     pace all pace coupling_coefficients.yace 1 0")

            # run

            lmp.command(f"thermo         100")
            #lmp.command(f"run {nsteps}")
            lmp.command(f"run 0")


        # declare compute pace variables

        dgradflag = 0
        run_lammps(dgradflag)
        lmp_pace = lmp.numpy.extract_compute("pace", LMP_STYLE_GLOBAL, LMP_TYPE_ARRAY)
        descriptor_grads = lmp_pace[ : len(atoms), : -1]
        return descriptor_grads
    #start = 'supercell_target.cif'
    file_prefix = 'iter_%d' % 0
    if type(start) == 'str':
        try:
            atoms = read(start)
        except:
            raise TypeError("unrecognized file type %s" % inp)
    elif type(start) == Atoms:
    #    except
        atoms = start 


    #atoms = read(start)

    write('%s.data' % file_prefix,atoms,format='lammps-data')
    start_arr = run_struct(atoms, '%s.data'% file_prefix)
    avg_start = vnp.average(start_arr,axis=0)
    var_start = vnp.var(start_arr,axis=0)
    vnp.save('target_descriptors.npy',avg_start)
    vnp.save('target_var_descriptors.npy',var_start)
    return avg_start, var_start
