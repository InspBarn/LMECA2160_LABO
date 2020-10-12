# *-* coding: utf-8 *-*

import math as mt
import numpy as np

import classCombustion as cbt

class MixedFuels:
	def __init__(self, fuel, gas, air_to_fuel_ratio):
		self.fuel = []
		self.fraction = np.zeros((len(fuel),))

		for (index,element) in fuel.items():
			(fuel_type, fuel_fraction) = element

			fuel_decomposition = {'O':0,'H':0,'C':0}
			for chem in fuel.index():
				i = fuel_type.index(chem) +1
				comp = 0
				while i < len(fuel_type) or i == -1:
					try:
						comp = comp*10 + int(fuel_type[i])
						i += 1
					except ValueError:
						i = -1
				fuel_decomposition[chem] = comp

			self.fuel = self.fuel.append(cbt.Combustion(
											x=fuel_decomposition{'O'},
											y=fuel_decomposition{'H'},
											z=fuel_decomposition{'C'}
										))
			self.fraction[index] = fuel_fraction

			self.reactants = fuel
			
			if 'N2' in self.fuel:
				self.reactants['N2'] = self.fuel['N2'] + air_to_fuel_ratio*(3.76*28)/(32+3.76*28)
			else:
				self.reactants['N2'] = air_to_fuel_ratio*(3.76*28)/(32+3.76*28)
			
			if 'O2' in self.fuel:
				self.reactants['O2'] = self.fuel['O2'] + air_to_fuel_ratio*32/(32+3.76*28)
			else:
				self.reactants['O2'] = air_to_fuel_ratio*32/(32+3.76*28)

			glob = sum(np.array(self.reactants.values()))
			for (chem,comp) in self.reactants:
				self.reactants[chem] = comp / glob

			self.products = gas
