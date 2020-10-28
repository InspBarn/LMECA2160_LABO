# *-* coding: utf-8 *-*

import math as mt
import numpy as np

import classCombustible as cbt

MH = 1.007
MC = 12.0107
MO = 15.9994
MN = 14.0067
rho_a = (.21e-3*MO*2 + .79e-3*MN*2) / (8.3144621*298.15/1.01325e5)

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
			fuel_decomposition = {'N':0,'O':0,'H':0,'C':0}
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
				w=fuel_decomposition['N'],
				x=fuel_decomposition['O'],
				y=fuel_decomposition['H'],
				z=fuel_decomposition['C']
			))
			self.fuel_fraction.append(ff)

		self.reactant = fuel

		if 'gas' in kwargs:
			self.product = kwargs["gas"]
		else:
			self.product = {}

		if 'AF_ratio' in kwargs:
			self.AF_ratio = kwargs["AF_ratio"]
		elif 'AFV_ratio' in kwargs:
			self.AF_ratio = kwargs["AFV_ratio"]  * rho_a/self.rho()
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
			w += ftype.w * ffraction
			x += ftype.x * ffraction
			y += ftype.y * ffraction
			z += ftype.z * ffraction
		return (w,x,y,z)

	# Return the molar mass in [kg/mol]
	def molar_mass(self):
		result = 0.0
		for (ftype, ffraction) in zip(self.fuel_type, self.fuel_fraction):
			result += ftype.molar_mass() * ffraction
		return result

	# Compute the mass fraction of fuel components and create a new attribute
	def _cpt_Mfraction(self):
		result = []
		for (ftype, ffraction) in zip(self.fuel_type, self.fuel_fraction):
			result.append(ffraction * ftype.molar_mass()/self.molar_mass())
		self.fuel_fraction_m = result
		return result

	# Compute the volumetric fraction of fuel components and create a new attribute
	def _cpt_Vfraction(self):
		result = []
		for (ftype, ffraction) in zip(self.fuel_type, self.fuel_fraction):
			result.append(ffraction * ftype.molar_mass()/self.molar_mass() / (ftype.rho()/self.rho()))
		self.fuel_fraction_v = result
		return result

	"""
	Return the mass density of the fuel [kg/m³]
	"""
	def rho(self):
		result = 0.0
		for (ftype, ffraction) in zip(self.fuel_type, self.fuel_fraction):
			result += ftype.rho() * ffraction
		return result

	"""
	Return the fuel LHV [J/mol]
	"""
	def LHV(self):
		result = 0.0
		for (ftype, ffraction) in zip(self.fuel_type, self.fuel_fraction):
			result += ftype.LHV() * ffraction
		return result

	"""
	Return the fuel LHV [J/kg]
	"""
	def LHV_m(self):
		return self.LHV() / self.molar_mass()

	"""
	Return the air to fuel ratio at stoechiometry
	"""
	def AF_stoech(self):
		if not hasattr(self, 'fuel_fraction_m'):
			self._cpt_Mfraction()
		result = 0.0
		for (ftype, ffraction) in zip(self.fuel_type, self.fuel_fraction_m):
			result += ftype.AF_stoech() * ffraction
		return result

	"""
	Return the volumetric air to fuel ratio at stoechiometry
	"""
	def AFV_stoech(self):
		if not hasattr(self, 'fuel_fraction_v'):
			self._cpt_Vfraction()
		result = 0.0
		for (ftype, ffraction) in zip(self.fuel_type, self.fuel_fraction_v):
			result += ftype.AFV_stoech() * ffraction
		return result

	"""
	Return the product to fuel ratio at stoechiometry
	"""
	def PF_stoech(self):
		if not hasattr(self, 'fuel_fraction_m'):
			self._cpt_Mfraction()
		result = 0.0
		for (ftype, ffraction) in zip(self.fuel_type, self.fuel_fraction_m):
			result += ftype.PF_stoech() * ffraction
		return result

	"""
	Return the volumetric product to fuel ratio at stoechiometry
	"""
	def PFV_stoech(self):
		if not hasattr(self, 'fuel_fraction_v'):
			self._cpt_Vfraction()
		result = 0.0
		for (ftype, ffraction) in zip(self.fuel_type, self.fuel_fraction_v):
			result += ftype.PFV_stoech() * ffraction
		return result

	"""
	Update the excess air coefficient of the combuxtion
	"""
	def excess_air_coefficient(self):
		self.lambda0 = self.AF_ratio / self.AF_stoech()

	def optimal_excess_air_coefficient(self):
		pass

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
