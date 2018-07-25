#==========================================================#
# Process: pp -> Z -> mu+mu-

# Author: Guojin Tseng
# Date: 2018.7.23
# Version: 1.0
#==========================================================#

"""
reference: https://matplotlib.org/xkcd/examples/api/histogram_path_demo.html

This example shows how to use a path patch to draw a bunch of
rectangles.  The technique of using lots of Rectangle instances, or
the faster method of using PolyCollections, were implemented before we
had proper paths with moveto/lineto, closepoly etc in mpl.  Now that
we have them, we can draw collections of regularly shaped objects with
homogeous properties more efficiently with a PathCollection.  This
example makes a histogram -- its more work to set up the vertex arrays
at the outset, but it should be much faster for large numbers of
objects
"""

#import guojin's module
import plot_card

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as path

class Gen_Histogram(object):

	def histogram(self, array, name):

		fig, ax = plt.subplots()
		n, bins = np.histogram(array, plot_card.bins)

		# get the corners of the rectangles for the histogram
		left = np.array(bins[:-1])
		right = np.array(bins[1:])
		bottom = np.zeros(len(left))
		top = bottom + n

		# we need a (numrects x numsides x 2) numpy array for the path helper
		# function to build a compound path
		XY = np.array([[left,left,right,right], [bottom,top,top,bottom]]).T

		# get the Path object
		barpath = path.Path.make_compound_path_from_polys(XY)

		# make a patch out of it
		patch = patches.PathPatch(barpath, facecolor='blue', edgecolor='gray', alpha=0.8)
		ax.add_patch(patch)

		# update the view limits
		ax.set_xlim(left[0], right[-1])
		ax.set_ylim(bottom.min(), top.max())
		plt.savefig("/home/guojintseng/Desktop/MC_Event_Generator_with_PDFs/fig/" + name + '.pdf')
