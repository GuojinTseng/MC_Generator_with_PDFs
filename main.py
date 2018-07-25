#==========================================================#
# Process: pp -> Z -> mu+mu-

# Author: Guojin Tseng
# Date: 2018.7.23
# Version: 1.0
#==========================================================#


#import guojin's module
import gen_histogram, event_generator, calc_xsection

class Generator(object):

	def __init__(self):
		self.xsection = calc_xsection.Calc_Xsection()
		self.events = event_generator.Event_Generator()
		self.histogram = gen_histogram.Gen_Histogram()

	def generate(self):
		weight_max = self.xsection.xsec()
		costh_list, SQ, Sy = self.events.generator(weight_max)
		self.histogram.histogram(SQ, 'Q')
		self.histogram.histogram(Sy, 'y')
		self.histogram.histogram(costh_list, 'costh')

if __name__ == "__main__":
	event = Generator().generate()
