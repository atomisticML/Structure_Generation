# creates:  nice.png

import numpy as np
import os
from ase import Atoms
from ase.build import bulk
from ase.io import write,read
from ase.utils import hsv

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from colorspacious import cspace_converter
from process_charges import *

def charge_arrays(ind,ats,typ='direct'):
    nats = len(ats)
    elems = [at.symbol for at in ats]
    chg_per_at = {symbol:[] for symbol in list(set(elems))}
    fix_chg_per_at = {symbol:[] for symbol in list(set(elems))}
    q_qeq,q_exact = get_charges(nats)
    q_qeq = -q_qeq.T
    print (np.shape(q_qeq))
    for ii in range(len(q_qeq)):
        for i,sym in enumerate(elems):
            chg_per_at[sym].append(q_qeq[ii][i])
    avg_per_at = {sym:np.average(chg_per_at[sym]) for sym in list(set(elems))}
    cmap = mpl.cm.get_cmap('PuOr_r')
    #cmap = mpl.cm.get_cmap('bwr')
    #cmap = mpl.cm.get_cmap('Spectral')
    if typ == 'direct':
        norm = mpl.colors.Normalize(vmin=np.min(q_qeq), vmax=np.max(q_qeq))
        return_q= q_qeq[ind]
        rgbas = [cmap(norm(qi)) for qi in return_q]
    elif typ == 'compare':
        q_diff = q_qeq
        all_reffed = []
        for qeqi in q_diff:
            all_reffed.append(qeqi)
        all_reffed = np.array(all_reffed)
        #q_diff = q_qeq[ind] -qref
        norm = mpl.colors.Normalize(vmin=np.min(all_reffed), vmax=np.max(all_reffed))
        #norm = mpl.colors.Normalize(vmin=-0.5, vmax=0.4)
        #return_q = q_diff[ind]
        return_q = all_reffed[ind]
        rgbas = [cmap(norm(qi)) for qi in return_q]
    elif typ == 'cold':
        q_diff = q_qeq
        #norm = mpl.colors.Normalize(vmin=np.min(q_diff), vmax=np.max(q_diff))
        norm = mpl.colors.Normalize(vmin=-0.4, vmax=0.4)
        return_q = q_diff[ind]
        rgbas = [cmap(norm(qi)) for qi in return_q]
    elif typ == 'dev':
        qref = np.array([avg_per_at[symi] for symi in elems]) 
        all_reffed = []
        for qeqi in q_qeq:
            all_reffed.append(qeqi)
        all_reffed = np.array(all_reffed)
        norm = mpl.colors.Normalize(vmin=np.min(all_reffed ), vmax= np.max(all_reffed))
        #norm = mpl.colors.Normalize(vmin=np.min(-0.4 ), vmax= 0.4)
        #norm = mpl.colors.Normalize(vmin=np.min(q_qeq - np.average(q_qeq,axis=0)), vmax=np.max(q_qeq - np.average(q_qeq,axis=0)))
        #norm = mpl.colors.Normalize(vmin=np.min(q_qeq[ind] - np.average(q_qeq[ind],axis=0)), vmax=np.max(q_qeq[ind] - np.average(q_qeq[ind],axis=0)))
        #return_q = q_qeq[ind] - np.average(q_qeq,axis=0)
        return_q = all_reffed[ind]
        print ('tot QEq', np.sum(q_qeq[ind]))
        #rgbas = [cmap(norm(qi)) for qi in return_q]

        rgbas_tmp = [cmap(norm(qi)) for qi in return_q]
        print ('Using amplified colors')
        rgbas = [(rgba[0]*3,rgba[1]*3,rgba[2]*3,rgba[3]) for rgba in rgbas_tmp]
    #if not os.path.isfile('./colorbar_%d.png'%ind):
    if not os.path.isfile('./colorbar.png'):
        fig,ax = plt.subplots()
        plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap))
        #plt.savefig('colorbar_%d.png' % ind,transparent=True)
        plt.savefig('colorbar.png',transparent=True)
    return return_q,rgbas

def plot_ats(f,ind = None,charge_color=False,typ='compare'):
    atoms = read(f)
    #new_cell = atoms.get_cell() * 1.2
    #atoms.set_cell(new_cell)
    #atoms.center()
    if ind == None:
        ind = int(f.split('.')[1])
    #typ = 'dev'
    #typ = 'compare'
    #typ = 'direct'
    prefix = '%s_%06d_ats' % (typ,ind)
    # view with ASE-GUI
    # view(atoms)
    #rotation = '10.5x, -12y, -1z'  # found using ASE-GUI menu 'view -> rotate'
    rotation = '41x, 76y, 40z'  # found using ASE-GUI menu 'view -> rotate'

    # Make colors
    color_map = {'U': [192/256,192/256,192/256],
            'Ni':[192/256,192/256,192/256],
            'Cr':[192/256,192/256,230/256],
            'O': [480/256, 1/256, 10/256],
            }
    tex_map =  {'U': 'ase2',
            'Ni':'ase2',
            'Cr':'ase2',
            'O':'glass',
            }
    if not charge_color:
        colors = [ color_map[atom.symbol] for atom in atoms]
    elif charge_color:
        charges,colors = charge_arrays(ind,atoms,typ=typ)
    tex = [tex_map[atom.symbol] for atom in atoms]

    # Textures
    #tex = ['jmol', ] * 288 + ['glass', ] * 288 + ['ase3', ] * 288 + ['vmd', ] * 288


    # Keywords that exist for eps, png, and povs
    generic_projection_settings = {
        'rotation': rotation,
        'colors': colors,
        'radii': None,
    }

    povray_settings = {  # For povray files only
        'display': False,  # Display while rendering
        'pause': False,  # Pause when done rendering (only if display)
        'transparent': True,  # Transparent background
        'canvas_width': 800,  # Width of canvas in pixels
        'canvas_height': None,  # Height of canvas in pixels
        'camera_dist': 50.,   # Distance from camera to front atom
        'image_plane': None,  # Distance from front atom to image plane
        # (focal depth for perspective)
        'camera_type': 'perspective',  # perspective, ultra_wide_angle
        'point_lights': [],             # [[loc1, color1], [loc2, color2],...]
        'area_light': [(2., 3., 40.),  # location
                       'White',       # color
                       .7, .7, 3, 3],  # width, height, Nlamps_x, Nlamps_y
        'background': 'White',        # color
        'textures': tex,  # Length of atoms list of texture names
        'celllinewidth': 0.05,  # Radius of the cylinders representing the cell
    }

    # Make flat png file
    #write('flat.png', atoms, **kwargs)

    # Make the color of the glass beads semi-transparent
    generic_projection_settings['colors'] = colors

    # Make the raytraced image
    # first write the configuration files, then call the external povray executable
    renderer = write('%s.pov'%prefix, atoms,
                     **generic_projection_settings,
                     povray_settings=povray_settings)
    renderer.render()

import glob


for ii,ind in enumerate(range(51)):
    f = 'sample.%d.cif' % ind
    #plot_ats('ats.10.cfg',charge_color=True)
    if ind == None:
        ind = int(f.split('.')[1])
    #typ = 'dev'
    typ = 'basic'
    #typ = 'direct'
    prefix = '%s_%06d_ats' % (typ,ind)
    if not os.path.isfile('%s.png' % prefix):
        plot_ats(f,charge_color=False,typ=typ)
