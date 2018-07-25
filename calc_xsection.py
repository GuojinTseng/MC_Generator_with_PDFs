#==========================================================#
# Process: pp -> Z -> mu+mu-

# Author: Guojin Tseng
# Date: 2018.7.23
# Version: 1.0
#==========================================================#

#import guojin's module
import calc_dsigma, param_card, run_card

from numba import jit

import math, os, sys
import random

import lhapdf
## initializes PDF member object (for protons)
lhapdf.initPDFSetByName("cteq6ll.LHpdf")
lhapdf.initPDF(0)

class Calc_Xsection(object):

	def __init__(self):
		self.dsig = calc_dsigma.Calc_Dsigma()


	def weight(self, hats, mu, x1, x2, costh_ii):
		# 1 for down-quarks, 2 for up, 3 for strange, 4 for charm and negative values for the corresponding anti-quarks. gluon is given by 21
		# up-type quarks
		qtype = 0
		w_ii = self.dsig.dsigma(costh_ii, hats, qtype) * ( lhapdf.xfx(x1, mu, 2) *  lhapdf.xfx(-x2, mu, 2) +  lhapdf.xfx(x1, mu, 4) *  lhapdf.xfx(-x2, mu, 4)  )
		w_ii = w_ii + self.dsig.dsigma(-costh_ii, hats, qtype) *  ( lhapdf.xfx(-x1, mu, 2) *  lhapdf.xfx(x2, mu, 2) +  lhapdf.xfx(-x1, mu, 4) *  lhapdf.xfx(x2, mu, 4) )
		# down-type quarks
		qtype = 1
		w_ii = w_ii + self.dsig.dsigma(costh_ii, hats, qtype) * ( lhapdf.xfx(x1, mu, 1) *  lhapdf.xfx(-x2, mu, 1) +  lhapdf.xfx(x1, mu, 3) *  lhapdf.xfx(-x2, mu, 3) )
		w_ii = w_ii + self.dsig.dsigma(-costh_ii, hats, qtype) * ( lhapdf.xfx(-x1, mu, 1) *  lhapdf.xfx(x2, mu, 1) +  lhapdf.xfx(-x1, mu, 3) *  lhapdf.xfx(x2, mu, 3) )
		# multiply by ranges and Jacobian, divide by the x1 and x2 since xfxaQ gives x * f(x)
		return w_ii

	# @jit
	def xsec(self):
		print '\n'
		print '----====================================================----'
		print 'pp --> Z --> mu+ mu-'
		print '----====================================================----'
		print '\n'

		random.seed(run_card.seed)

		print "hadron com energy:", run_card.ECM, "GeV"

		# choose the transform mass and width
		MTR = param_card.Q_min
		GammaTR = param_card.Q_min

		# we also need the "ranges" (i.e. x2 - x1)
		# for costh this is 1 - (-1) = 2
		deltath = 2
		# choose rho limits (see transformation)
		rho1 = math.atan( param_card.Q_min**2 - MTR**2 ) / (GammaTR*MTR)
		rho2 = math.atan( (run_card.s - MTR**2) / (GammaTR*MTR) )
		deltarho = rho2 - rho1


		# loop over N phase space points,
		# sum the weights up (sum_w) and the weights squared (sum_w_sq) (for the variance)
		sum_w = 0
		sum_w_sq = 0

		# also define maximum point variable
		w_max = 0
		costh_max = 0
		Q_max = 0

		print 'integrating for cross section and maximum!'
		print '...'
		for ii in range(0, run_card.N):
			# random costheta and rho
			costh_ii = -1 + random.random() * deltath
			rho = rho1 + random.random() * deltarho
			# jacobian
			Jac = (MTR * GammaTR) / ( math.cos(rho)**2 * run_card.s)
			# get s hat
			hats = MTR * GammaTR * math.tan(rho) + MTR**2
			# get maximum rapidity of dilepton, Y and find the range of integration for y (=2*Y)
			Y = - 0.5 * math.log(hats/run_card.s)
			deltay = 2 * Y
			# get a random value of y
			y = ( (2 * random.random()) - 1 ) * Y
			# calculate momentum fractions x1, x2
			x1 = math.sqrt(hats/run_card.s) * math.exp(y)
			x2 = math.sqrt(hats/run_card.s) * math.exp(-y)
			# get the scale Q
			Q = math.sqrt(hats)
			# set the scale for the pdfs
			mu = param_card.MZ
			# calc. phase space point weight
			w_ii = self.weight(hats, mu, x1, x2, costh_ii) * deltath * deltarho * deltay * Jac / (x1 * x2)
			# add to the sums
			sum_w = sum_w + w_ii
			sum_w_sq = sum_w_sq + w_ii**2
			# check if higher than maximum
			if w_ii > w_max:
				w_max = w_ii
				costh_max = costh_ii
				Q_max = math.sqrt(hats)

		 # calculate cross section
		sigma = sum_w / run_card.N

		# and its error through the variance
		variance = sum_w_sq/run_card.N - (sum_w/run_card.N)**2
		error = math.sqrt(variance/run_card.N)

		print 'done integrating!'
		print '\n'
		print 'maximum value of dsigma = ', w_max, 'found at costh = ', costh_max
		print 'total cross section =', sigma * param_card.pb_convert, '+-', error * param_card.pb_convert, 'pb'

		return w_max
