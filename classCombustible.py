# *-* coding: utf-8 *-*

import math as mt
import numpy as np

MH = 1.007
MC = 12.0107
MO = 15.9994

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

	# Return the molar mass in [kg/mol]
	def molar_mass(self):
		return (self.x*MO + self.y*MH + self.z*MC) * 1e-3

	def LHV(self):
		if (self.x == self.y == self.z == 0) \
			or (self.z==1 && self.y==0 && self.x==2) \
			or (self.z==0 && self.y==2 && self.x==1):
			# No Fuel or CO2 or H2O
			return 0
		elif (self.y==2) && (self.x==self.z==0):
			# Gaseous Fuell // H2
			return 241800 / self.molar_mass()
		elif (self.z==1) && (self.y==4) && (self.x==0):
			# Gaseous Fuel // CH4
			return 802400 / self.molar_mass()
		elif (self.z==3) && (self.y==8) && (self.x==0):
			# Gaseous Fuel // C3H8
			return 2044400 / self.molar_mass()
		elif (self.z==self.x==1) && (self.y==0):
			# Gaseous Fuel // Carbon monoxyde
			return 282400 / self.molar_mass()
		elif (self.z == 2*self.y+2) && (self.x == 0):
			# Liquid Fuel // Paraffin
			return (612170*self.z - 207820) / self.molar_mass()
		elif (self.z == 2*self.y) && (self.x == 0):
			# Liquid Fuel // Naphtene
			return 612170*self.z / self.molar_mass()
		elif (self.z == 2*self.y-2) && (self.x == 0):
			# Liquid Fuel // Diolefin
			return (612170*self.z - 84184) / self.molar_mass()
		elif (self.z == self.y == 6) && (self.x == 0):
			# Liquid Fuel // Benzene
			return 6 * 522200 / self.molar_mass()
		elif self.x==0:
			# Solid Fuel // Coke
			return (393400 + 131760*self.y) / self.molar_mass()
		else:
			# Solid Fuel // General CHyOx
			return (393400 + 102550*self.y - (111000+102250*self.y)*self.x/(1+self.y/2.)) / self.molar_mass()


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
