# *-* coding: utf-8 *-*

import math as mt
import numpy as np

class Combustion:
	def __init__(x: int=0, y: int=4, z: int=1, LHV: float=50E6, AF: float=25.74):
		self.x = x
		self.y = y
		self.z = z

		self.LHV = LHV
		self.AF = AF

	def AF_stoech(self):
		return (32.0 + 3.76*28) * (z + (y-2*x)/4) / (12*x + y + 16*x)

	def AFV_stoech(self):
		return 4.76 * (z + (y-2*x)/4)

	def PF_stoech(self):
		#return (44*z + 9*y + 3.76*28*(z+(y-2*x)/4)) / (12*z + y + 16*x)
		return self.AF_stoech() + 1

	def PFV_stoech(self):
		return z + y/2 + 3.76*(z + (y-2*x)/4)

	def carbon_monoxide_coefficient(self):
		return self.CO / (self.CO + self.CO2)

	def excess_air_coefficient(self):
		self.lambda0 = self.AF / self.AF_stoech()

	def equivalence_ratio(self):
		return self.lambda0**(-1)

	def balance_equation_to_fuel(eps: float=1E-3):
		left_hand_side = 4.76*self.O2 + (4.76 + 3.76*(y-2*x)/(4*z))*self.CO2 + (2.88 + 3.76*(y-2*x)/(4*z) - 0.88*y/(4*z))*self.CO
		return (left_hand_side - 1 < eps)
