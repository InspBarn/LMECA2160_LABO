# *-* coding: utf-8 *-*

import math as mt
import numpy as np

MH = 1.007
MC = 12.0107
MO = 15.9994
MN = 14.0067

hf_fuel = [x*1e3 for x in [-74.8, -84.7, -103.8, -126.2]] # [J/mol]
hf_CO2  = -393.52e3
hf_H20  = -241.8e3

"""
Combustible est un objet représentant un combustible ainsi que sa combustion. Un combustible
est defini par la formule :
                                          Cz Hy Ox
Les paramètres sont :
	- w : indice d'azote
	- x : indice d'oxygene
	- y : indice d'hydrogene
	- z : indice de carbone

Le méthane (CH4) est le combustible defini par defaut.
"""
class Combustible:
	def __init__(self, w: int=0, x: int=0, y: int=4, z: int=1):
		self.w = w
		self.x = x
		self.y = y
		self.z = z

	# Return the molar mass in [kg/mol]
	def molar_mass(self):
		return (self.w*MN + self.x*MO + self.y*MH + self.z*MC) * 1e-3

	# We assume we have a gas, here (clearly do not want to do more researches !)
	def rho(self):
		return self.molar_mass() / (8.3144621*298.15/1.01325e5)

	"""
	Return the Lower Heating Value (LHV) of the fuel [J/mol]
	"""
	def LHV(self):
		if (self.x == self.y == self.z == 0) \
			or (self.z==1 and self.y==0 and self.x==2) \
			or (self.z==0 and self.y==2 and self.x==1):
			# No Fuel or CO2 or H2O
			return 0
		elif (self.y==2) and (self.x==self.z==0):
			# Gaseous Fuell // H2
			return 241800
		elif (self.z==self.x==1) and (self.y==0):
			# Gaseous Fuel // Carbon monoxyde
			return 282400
		else:
			# Gaseous Fuel // Carbon oxydes CnHm
			return hf_fuel[self.z-1] - self.z*hf_CO2 - self.y/2.*hf_H20

	"""
	Return the air to fuel ratio of the combustible at stoechiometry
	"""
	def AF_stoech(self):
		if (self.x==self.y==self.z==0):
			return 0
		return (2*MO + 3.76*2*MN) * (self.z+(self.y-2*self.x)/4) / (self.z*MC + self.y*MH + self.x*MO)

	"""
	Return the volumetric air to fuel ratio of the combustible at stoechiometry
	"""
	def AFV_stoech(self):
		if (self.x==self.y==self.z==0):
			return 0
		return 4.76 * (self.z + (self.y-2*self.x)/4)

	"""
	Return the product to fuel ratio of the combustible at stoechimotry
	"""
	def PF_stoech(self):
		if (self.x==self.y==self.z==0):
			return 1
		return (self.z*(MC+2*MO) + self.y/2*(2*MH+MO) + 3.76*2*MN*(self.z+(self.y-2*self.x)/4)) / (self.z*MC + self.y*MH + self.x*MO)
		#return self.AF_stoech() + 1

	"""
	Return the volumetric product to fuel ratio of the combustible at stoechiometry
	"""
	def PFV_stoech(self):
		return self.z + self.y/2 + 3.76*(self.z + (self.y-2*self.x)/4)

"""
	def balance_equation_to_fuel(eps: float=1E-3):
		left_hand_side = 4.76*self.O2 + (4.76 + 3.76*(y-2*x)/(4*z))*self.CO2 + (2.88 + 3.76*(y-2*x)/(4*z) - 0.88*y/(4*z))*self.CO
		return (left_hand_side - 1 < eps)
"""
