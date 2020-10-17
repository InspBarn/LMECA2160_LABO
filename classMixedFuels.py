# *-* coding: utf-8 *-*

import math as mt
import numpy as np

import classCombustible as cbt

"""
MixedFuels est un objet représentant un fuel le plus général possible consistant en un
mélange de 'N' composants chimiques différents de fraction 'X' dont les paramètres sont :
	- fuel_type     : liste d'objets 'Combustible' représentant la combustion de chaque composants
						du fuel. Notons que chaque composants ne brule pas spcialement (e.g. 'N2')
	- fuel_fraction : np.array de floats representant la fraction de chaque composant. La somme
						des éléments vaut 1.
	- reactant      : dictionnaire des réactifs {"comp": frac} ou :
						* comp est le nom du composant
						* frac est la fraction massique du composants dans les réactifs
	- product       : dictionnaire des produits {"comp": frac} ou :
						* comp est le nom du composant
						* frac est la fraction massique du composants dans les produits
	- AF_ratio      : air to fuel ratio, float

Le fuel defini par defaut est le methane (CH4) sans melange. Le ratio d'air au fuel est defini
a la stoechiometrie.
"""
class MixedFuels:
	def __init__(self, fuel: dict={"CH4": 1.0}, **kwargs):
		self.fuel_type = []
		self.fuel_fraction = []
		for (ft,ff) in fuel.items():
			fuel_decomposition = {'O':0,'H':0,'C':0}
			for chem in fuel_decomposition.keys():
				if chem in ft:
					i,c = ft.index(chem) +1, 0
					while i != -1:
						try:
							c = c*10 + int(ft[i])
							i += 1
						except:
							c = 1 if c==0 else c
							i = -1
					fuel_decomposition[chem] = c

			self.fuel_type.append(cbt.Combustible(
				x=fuel_decomposition['O'],
				y=fuel_decomposition['H'],
				z=fuel_decomposition['C']
			))
			self.fuel_fraction.append(ff)

		self.reactant = fuel

		if kwargs.get("gas"):
			self.product  = kwargs["gas"]
		else:
			self.product = {}

		if kwargs.get("AF_ratio"):
			self.AF_ratio = kwargs["AF_ratio"]
		else:
			self.AF_ratio = self.AF_stoech()

		if 'N2' in fuel:
			self.reactant['N2'] = fuel['N2'] + self.AF_ratio*(3.76*28)/(32+3.76*28)
		else:
			self.reactant['N2'] = self.AF_ratio*(3.76*28)/(32+3.76*28)

		if 'O2' in fuel:
			self.reactant['O2'] = fuel['O2'] + self.AF_ratio*32/(32+3.76*28)
		else:
			self.reactant['O2'] = self.AF_ratio*32/(32+3.76*28)

		glob = sum(self.reactant.values())
		for (chem,comp) in self.reactant.items():
			self.reactant[chem] = comp / glob

		glob = sum(self.product.values())
		for (chem,comp) in self.product.items():
			self.product[chem] = comp / glob

		self.excess_air_coefficient()

	def get_xyz(self):
		x,y,z = 0,0,0
		for (ftype, ffraction) in zip(self.fuel_type, self.fuel_fraction):
			x += ftype.x * ffraction
			y += ftype.y * ffraction
			z += ftype.z * ffraction
		return (x,y,z)

	"""
	Return the mass density of the fuel [kg/m³]
	"""
	def rho(self):
		result = 0.0
		for (index,fuel_type) in enumerate(self.fuel_type):
			result += fuel_type.rho() * self.fuel_fraction[index]
		return result

	"""
	Return the fuel LHV [J/kg]
	"""
	def LHV(self):
		result = 0.0
		for (ftype, ffraction) in zip(self.fuel_type, self.fuel_fraction):
			result += ftype.LHV() * ffraction
		return result

	"""
	Return the air to fuel ratio at stoechiometry
	"""
	def AF_stoech(self):
		result = 0.0
		for (ftype, ffraction) in zip(self.fuel_type, self.fuel_fraction):
			result += ftype.AF_stoech() * ffraction
		return result

	"""
	Return the volumetric air to fuel ratio at stoechiometry
	"""
	def AFV_stoech(self):
		result = 0.0
		for (index,fuel_type) in enumerate(self.fuel_type):
			result += fuel_type.AFV_stoech() * self.fuel_fraction[index]
		return result

	"""
	Return the product to fuel ratio at stoechiometry
	"""
	def PF_stoech(self):
		result = 0.0
		for (index,fuel_type) in enumerate(self.fuel_type):
			result += fuel_type.PF_stoech() * self.fuel_fraction[index]
		return result

	"""
	Return the volumetric product to fuel ratio at stoechiometry
	"""
	def PFV_stoech(self):
		result = 0.0
		for (index,fuel_type) in enumerate(self.fuel_type):
			result += fuel_type.PFV_stoech() * self.fuel_fraction[index] 
		return result

	"""
	Update the excess air coefficient of the combuxtion
	"""
	def excess_air_coefficient(self):
		self.lambda0 = self.AF_ratio / self.AF_stoech()

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
