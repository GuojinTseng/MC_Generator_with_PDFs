#==========================================================#
# Process: pp -> Z -> mu+mu-

# Author: Guojin Tseng
# Date: 2018.7.23
# Version: 1.0
#==========================================================#

#import guojin's module
import param_card, run_card

import math

class Calc_Dsigma(object):

	def dsigma(self, costh, hats, qtype):
		# constants and functions that appear in the differential cross section
		kappa = math.sqrt(2) * param_card.G_Fermi * param_card.MZ**2 / (4 * math.pi * param_card.alpha)
		chi1 = kappa * hats * ( hats - param_card.MZ**2 ) / (  (hats-param_card.MZ**2)**2 + param_card.GAMMAZ**2*param_card.MZ**2 )
		chi2 = kappa**2 * hats**2 / (  (hats-param_card.MZ**2)**2 + param_card.GAMMAZ**2*param_card.MZ**2 )
		# CL and CR for leptons
		CVe = -0.5 + 2 * param_card.sw2
		CAe = -0.5

		# up-type quarks
		if qtype is 0:
			CVf = 0.5 - (4./3.) * param_card.sw2
			CAf = 0.5
			Qf = 2./3. # charge of the quark
		# down-type quarks
		if qtype is 1:
			CVf = - 0.5 + (2./3.) * param_card.sw2
			CAf = -0.5
			Qf = -1./3.
		A0 = Qf**2 - 2 * Qf * CVe * CVf * chi1 + (CAe**2 + CVe**2) * (CAf**2 + CVf**2) * chi2  # see notes
		A1 = -4 * Qf * CAe * CAf * chi1 + 8 * CAe * CVe * CAf * CVf * chi2
		PREFAC = (2. * math.pi) * param_card.alpha**2 / (4.*hats) /3. # 2 * pi comes from d-phi integral, 1/3 from initial colour averaging
		return  PREFAC * ( A0 * ( 1 + costh**2 ) + A1 * costh )
