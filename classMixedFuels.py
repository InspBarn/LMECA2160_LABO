# *-* coding: utf-8 *-*

import math as mt
import numpy as np

import classCombustible as cbt

"""
MixedFuels est un objet représentant un fuel le plus général possible consistant en un
mélange de 'N' composants chimiques différents de fraction 'X' dont les paramètres sont :
	- fuelT    : liste d'objets 'Combustible' représentant la combustion de chaque composants
				du fuel. Notons que chaque composants ne brule pas spcialement (e.g. 'N2')
	- fraction : np.array de floats representant la fraction de chaque composant. La somme
				des éléments vaut 1.
	- reactant : dictionnaire des réactifs {"comp": frac} ou :
				* comp est le nom du composant
				* frac est la fraction massique du composants dans les réactifs
	- product  : dictionnaire des produits {"comp": frac} ou :
				* comp est le nom du composant
				* frac est la fraction massique du composants dans les produits
	- AF_ratio : air to fuel ratio, float

Le fuel defini par defaut est le methane (CH4) sans melange. Le ratio d'air au fuel est defini
a la stoechimotrie.
"""
class MixedFuels:
	def __init__(self, fuel: dict, gas: dict, AF_ratio: float):
		fuelT, fraction = [],[]
		for (fuel_type,fuel_fraction) in fuel.items():
			fuel_decomposition = {'O':0,'H':0,'C':0}
			for chem in fuel.index():
				i,c = fuel_type.index(chem) +1, 0
				while (i<len(fuel_type)) or (i==-1):
					try:
						c = c*10 + int(fuel_type[i])
						i += 1
					except ValueError:
						i = -1
				fuel_decomposition[chem] = comp

			fuelT.append(cbt.Combustible(
				x=fuel_decomposition{'O'},
				y=fuel_decomposition{'H'},
				z=fuel_decomposition{'C'}
			))
			fraction.append(fuel_fraction)

		self.fuelT    = fuelT
		self.fraction = np.array(fraction)

		self.reactant = fuel
		self.product  = gas

		if 'N2' in self.fuel:
			self.reactant['N2'] = self.fuel['N2'] + AF_ratio*(3.76*28)/(32+3.76*28)
		else:
			self.reactant['N2'] = AF_ratio*(3.76*28)/(32+3.76*28)

		if 'O2' in self.fuel:
			self.reactant['O2'] = self.fuel['O2'] + AF_ratio*32/(32+3.76*28)
		else:
			self.reactant['O2'] = AF_ratio*32/(32+3.76*28)

		glob = sum(self.reactant.values())
		for (chem,comp) in self.reactant.items():
			self.reactant[chem] = comp / glob

		glob = sum(self.product.values())
		for (chem,comp) in self.product.items():
			self.product[chem] = comp / glob

		self.AF_ratio = AF_ratio
		self.lambda0  = None

		return self

	"""
	Return the fuel LHV
	"""
	def LHV(self):
		result = 0.0
		for (index,fuel_type) in enumerate(self.fuelT):
			result += fuel_type.LHV() * self.fraction[index]
		return result

	"""
	Return the air to fuel ratio at stoechiometry
	"""
	def AF_stoech(self):
		result = 0.0
		for (index,fuel_type) in enumerate(self.fuelT):
			result += fuel_type.AF_stoech() * self.fraction[index]
		return result

	"""
	Return the volumetric air to fuel ratio at stoechiometry
	"""
	def AFV_stoech(self):
		result = 0.0
		for (index,fuel_type) in enumerate(self.fuelT):
			result += fuel_type.AFV_stoech() * self.fraction[index]
		return result

	"""
	Return the product to fuel ratio at stoechiometry
	"""
	def PF_stoech(self):
		result = 0.0
		for (index,fuel_type) in enumerate(self.fuelT):
			result += fuel_type.PF_stoech() * self.fraction[index]
		return result

	"""
	Return the volumetric product to fuel ratio at stoechiometry
	"""
	def PFV_stoech(self):
		result = 0.0
		for (index,fuel_type) in enumerate(self.fuelT):
			result += fuel_type.PFV_stoech() * self.fraction[index] 
		return result

	"""
	Update the excess air coefficient of the combuxtion
	"""
	def excess_air_coefficient(self):
		self.lambda0 = self.AF / self.AF_stoech()

	"""
	Update the CO coefficient
	"""
	def carbon_monoxide_coefficient(self):
		self.k = self.product['CO'] / (self.product['CO'] + self.product['CO2'])

	"""
	Return the equivalence ratio
	"""
	def equivalence_ratio(self):
		try:
			return self.lambda0**(-1)
		except ValueError:
			self.excess_air_coefficient()
			return self.lambda0**(-1)
