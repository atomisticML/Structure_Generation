import numpy as np
import math
import matplotlib.pyplot as plt
from gen_labels import *


def alter_axis(ax,x,y):
	minx = np.amin(x)
	maxx = np.amax(x)
	miny = 0.0 #np.amin(y)
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


ntemps =3
hist_range = range(ntemps)
combined_hist = np.array([])

for ih in hist_range:
	this_hist = np.load('full_avg_hist_%d.npy'%ih)
	if np.shape(combined_hist) == (0,):
		combined_hist = this_hist
	elif np.shape(combined_hist) != (0,):
		combined_hist = np.append(combined_hist,this_hist,axis=0)

tol = 2

combined_hist = combined_hist.T

for idescrow in range(np.shape(combined_hist)[0]):
	this_combined_hist= combined_hist[idescrow]

	this_combined_hist = np.around(this_combined_hist,tol)

	#binned = np.bincount(combined_hist)
	#binned,edges = np.histogram(this_combined_hist,range=(min(this_combined_hist),max(this_combined_hist)), bins = 11)
	binned,edges = np.histogram(this_combined_hist,range=(min(flatten(combined_hist)),max(flatten(combined_hist))), bins = 11)


	fig,ax = plt.subplots()

	xax = edges[1:]

	xax = xax -1

	print (sorted(this_combined_hist))

	print (binned)
	print (edges)

	#plt.bar( edges ,binned,color='None',edgecolor='b',clip_on = False)
	plt.bar( xax, binned, width=0.5, color='None',edgecolor='b',clip_on = False)

	plt.xlabel("Value",fontsize=14)
	plt.ylabel(r"Frequency",fontsize=14)

	alter_axis(ax,xax,binned)
	plt.tight_layout()
	plt.savefig('histogram_%d.png' % idescrow,transparent=True)

	fig.clf()
