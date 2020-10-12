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

			fuelT.append(cbt.Combustion(
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
		for (chem,comp) in self.reactant:
			self.reactant[chem] = comp / glob

		return self
