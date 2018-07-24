#==========================================================#
# Process: pp -> Z -> mu+mu-

# Author: Guojin Tseng
# Date: 2018.7.23
# Version: 1.0
#==========================================================#

#import guojin's module
import run_card, param_card

class Event_Output(object):

	def output_headers(self):
		with open("event.txt", "a") as events:
			events.write("#**************************************************")
			events.write("\n")
			events.write("#*            Guojin's Event Generator            *")
			events.write("\n")
			events.write("#*                                                *")
			events.write("\n")
			events.write("#*           *                       *            *")
			events.write("\n")
			events.write("#*             *                   *              *")
			events.write("\n")
			events.write("#*               * * G-J Tseng * *                *")
			events.write("\n")
			events.write("#*             *                   *              *")
			events.write("\n")
			events.write("#*           *                       *            *")
			events.write("\n")
			events.write("#*                                                *")
			events.write("\n")
			events.write("#*                                                *")
			events.write("\n")
			events.write("#*    VERSION 1.0                 2018-07-23      *")
			events.write("\n")
			events.write("#*                                                *")
			events.write("\n")
			events.write("#*                  Find me at                    *")
			events.write("\n")
			events.write("#*             guojintseng@pku.edu.cn             *")
			events.write("\n")
			events.write("#*                                                *")
			events.write("\n")
			events.write("#**************************************************")
			events.write("\n")
			events.write("#*                                                *")
			events.write("\n")
			events.write("#*              All rights reserved.              *")
			events.write("\n")
			events.write("#*                                                *")
			events.write("\n")
			events.write("#**************************************************")
			events.write("\n")
			events.write("\n")
			events.write("#*********************run_card*********************")
			events.write("COM Energy in GeV: ", run_card.ECM)
			events.write("Random Seed: ", run_card.seed )
			events.write("Number of Sprinkled Points: ", run_card.N)
			events.write("Number of Generated events: ", run_card.Nevents )
			events.write("\n")
			events.write("\n")
			events.write("#*******************param_card*********************")
			events.write("Z Boson Mass: ", param_card.MZ)
			events.write("Z Boson Width: ", param_card.GAMMAZ)
			events.write("ALPHA_EM @MZ: ", param_card.alpha)
			events.write("G_Fermi Constant: ", param_card.G_Fermi)
			events.write("# sin^2(weinberg angle)", param_card.sw2)

	def output(self, i, p1, p2, p3, p4):
		with open("event.txt","a") as events:
			events.write("===================="+"event "+str(i)+"====================")
			events.write("\n")
			events.write("pq1:	"+str(p1))
			events.write("\n")
			events.write("pq1:	"+str(p2))
			events.write("\n")
			events.write("pmm:	"+str(p3))
			events.write("\n")
			events.write("pmp:	"+str(p4))
			events.write("\n")
			events.write("\n")
