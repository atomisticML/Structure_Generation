import numpy as np
import glob
import json
from ase.io import read,write
import math

def alter_axis(ax,x,y):
	minx = 0 #np.amin(x)
	maxx = np.amax(x)
	miny =  0.0 #np.amin(y)
	maxy = np.amax(y)
	x_inc = 4. #0.0005
	y_inc = 0.05#*1000
	plt.axis((math.floor(minx/x_inc)/(1/x_inc),math.ceil(maxx/x_inc)/(1/x_inc),
		math.floor(miny/y_inc)/(1/y_inc),math.ceil(maxy/y_inc)/(1/y_inc)))
	ax.spines['bottom'].set_position(('outward',2))
	ax.spines['left'].set_position(('outward',10))
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	return

import matplotlib.pyplot as plt

fig,ax = plt.subplots()


all_resids=[]
separated_resids = []
for iresid in sorted(glob.glob('residuals*.npy')):
	tmpresid = np.load(iresid)
	separated_resids.append(tmpresid)
	tmpresid = list(tmpresid)
	all_resids.extend(tmpresid)
x_ax = list(range(len(all_resids)))

#used_inds = []
#for sr in separated_resids:
#	plt.scatter(x_ax

plt.scatter(range(len(all_resids)),all_resids,color='None',edgecolor='b',clip_on = False)

plt.xlabel("MC Step",fontsize=14)
plt.ylabel(r"Summed descriptor residuals",fontsize=14)

alter_axis(ax,range(len(all_resids)),all_resids)
plt.tight_layout()
plt.savefig('combined_residuals.png',transparent=True)

fig.clf()

