#==========================================================#
# Process: pp -> Z -> mu+mu-

# Author: Guojin Tseng
# Date: 2018.7.23
# Version: 1.0
#==========================================================#

#import guojin's module
import calc_dsigma, event_output, run_card, param_card

import math, random, os, sys

import lhapdf
## initializes PDF member object
lhapdf.initPDFSetByName("cteq6ll.LHpdf")
lhapdf.initPDF(0)


class Event_Generator(object):

	def __init__(self):
		self.dsig = calc_dsigma.Calc_Dsigma()
		self.eoutput = event_output.Event_Output()


	def boost(self, fourvector, beta):
		gamma = math.sqrt(1 / (1 - beta**2))
		boosted_vector = [gamma * fourvector[0] - gamma * beta * fourvector[3], fourvector[1], fourvector[2], - gamma * beta * fourvector[0] + gamma * fourvector[3]]
		return boosted_vector


	def weight(self, hats, mu, x1, x2, costh_i):
		'''
		1 for down-quarks, 2 for up, 3 for strange,
		4 for charm and negative values for the corresponding anti-quarks. gluon is given by 21
		'''
		# up-type quarks
		qtype = 0
		w_i = self.dsig.dsigma(costh_i, hats, qtype) * ( lhapdf.xfx(x1, mu, 2) *  lhapdf.xfx(-x2, mu, 2) +  lhapdf.xfx(x1, mu, 4) *  lhapdf.xfx(-x2, mu, 4)  )
		w_i = w_i + self.dsig.dsigma(-costh_i, hats, qtype) *  ( lhapdf.xfx(-x1, mu, 2) *  lhapdf.xfx(x2, mu, 2) +  lhapdf.xfx(-x1, mu, 4) *  lhapdf.xfx(x2, mu, 4) )
		# down-type quarks
		qtype = 1
		w_i = w_i + self.dsig.dsigma(costh_i, hats, qtype) * ( lhapdf.xfx(x1, mu, 1) *  lhapdf.xfx(-x2, mu, 1) +  lhapdf.xfx(x1, mu, 3) *  lhapdf.xfx(-x2, mu, 3) )
		w_i = w_i + self.dsig.dsigma(-costh_i, hats, qtype) * ( lhapdf.xfx(-x1, mu, 1) *  lhapdf.xfx(x2, mu, 1) +  lhapdf.xfx(-x1, mu, 3) *  lhapdf.xfx(x2, mu, 3) )
		return w_i


	def generator(self, w_max):
		delta_th = 2
		MTR = param_card.Q_min
		GammaTR = param_card.Q_min

		# choose rho limits
		rho1 = math.atan( param_card.Q_min**2 - MTR**2 ) / (GammaTR*MTR)
		rho2 = math.atan( (run_card.s - MTR**2) / (GammaTR*MTR) )
		delta_rho = rho2 - rho1
		random.seed(run_card.seed)

		# Generating unweighted events according to the 'hit-or-miss' method. Maybe there are another methods!
		print 'generating events...'

		costh_list = []
		SQ = []
		Sy = []

		# counter of events generate
		i = 0

		self.eoutput.output_headers()

		while i < run_card.Nevents:
			# generate random cos(theta) and rho
			costh_i = -1 + random.random() * delta_th
			rho = rho1 + random.random() * delta_rho

			# calculate jacobian factor
			Jacobian_Factor= (MTR * GammaTR) / ( math.cos(rho)**2 * run_card.s)
			# calculate hats
			hats = MTR * GammaTR * math.tan(rho) + MTR**2
			# get maximum rapidity of dilepton, Y and find the range of integration for y (=2*Y)
			Y = - 0.5 * math.log(hats/run_card.s)
			delta_y = 2 * Y
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
			w_i = self.weight(hats, mu, x1, x2, costh_i) * delta_th * delta_rho * delta_y * Jacobian_Factor/ (x1 * x2)
			# now divide by maximum and compare to pro bability
			prob = w_i / w_max
			rand_num = random.random()
			# if the random number is less than the probability of the S point
			# accept

			n = i/5000

			if rand_num < prob:
				i = i + 1

				costh_list.append(costh_i)
				SQ.append(Q)
				Sy.append(y)

				beta = (x2 - x1) / (x2 + x1)
				# generate random phi
				phi = random.random() * 2 * math.pi
				sinphi = math.sin(phi)
				cosphi = math.cos(phi)
				sinth = math.sqrt( 1 - costh_i**2 )
				pq1 = [ 0.5 * x1 * run_card.ECM, 0., 0., 0.5 * x1 * run_card.ECM ]
				pq2 = [ 0.5 * x2 * run_card.ECM, 0., 0., - 0.5 * x2 * run_card.ECM ]
				pem = [ 0.5 * Q, 0.5 * Q * sinth * cosphi, 0.5 * Q * sinth * sinphi, 0.5 * Q * costh_i ]
				pep = [ 0.5 * Q, - 0.5 * Q * sinth * cosphi, - 0.5 * Q * sinth * sinphi, - 0.5 * Q * costh_i ]
				# boost them to the lab frame
				pem_boosted = self.boost(pem,beta)
				pep_boosted = self.boost(pep,beta)

				# print 'q1', pq1
				# print 'q2', pq2
				# print 'mu-', pem_boosted
				# print 'mu+', pep_boosted
				# print '\n'

				self.eoutput.output(i, pq1, pq2, pem, pep)

				if n < i/5000:
					print 100. * float(i)/run_card.Nevents, "%% events has been generated."

		return costh_list, SQ, Sy
