# *-* coding: utf-8 *-*

import math as mt
import numpy as np

"""
Combustible est un objet représentant un combustible ainsi que sa combustion. Un combustible
est defini par la formule :
                                          Cz Hy Ox
Les paramètres sont :
	- x : indice d'oxygene
	- y : indice d'hydrogene
	- z : indice de carbone

Le méthane (CH4) est le combustible defini par defaut.
"""
class Combustible:
	def __init__(x: int=0, y: int=4, z: int=1):
		self.x = x
		self.y = y
		self.z = z

	"""
	Return the air to fuel ratio of the combustible at stoechiometry
	"""
	def AF_stoech(self):
		return (32.0 + 3.76*28) * (z + (y-2*x)/4) / (12*x + y + 16*x)

	"""
	Return the volumetric air to fuel ratio of the combustible at stoechiometry
	"""
	def AFV_stoech(self):
		return 4.76 * (z + (y-2*x)/4)

	"""
	Return the product to fuel ratio of the combustible at stoechimotry
	"""
	def PF_stoech(self):
		#return (44*z + 9*y + 3.76*28*(z+(y-2*x)/4)) / (12*z + y + 16*x)
		return self.AF_stoech() + 1

	"""
	Return the volumetric product to fuel ratio of the combustible at stoechiometry
	"""
	def PFV_stoech(self):
		return z + y/2 + 3.76*(z + (y-2*x)/4)

"""
	def balance_equation_to_fuel(eps: float=1E-3):
		left_hand_side = 4.76*self.O2 + (4.76 + 3.76*(y-2*x)/(4*z))*self.CO2 + (2.88 + 3.76*(y-2*x)/(4*z) - 0.88*y/(4*z))*self.CO
		return (left_hand_side - 1 < eps)
"""
