# *-* coding: utf-8 *-*

import csv
import math as mt
import numpy as np
import pandas as pd

import classCombustible as cbt
import classMixedFuels as mxf

with open('Labo_Result_2020_09_28_14_21.txt') as labo_result:
	reader = csv.reader(labo_result)
	for i,row in enumerate(reader):
		row = row[0].split('\t')
		row = np.expand_dims([float(x) for x in row], axis=0)
		try:
			results = np.concatenate((results,row))
		except NameError:
			if i > 3:
				results = row

results = pd.DataFrame(results)
results.columns = ['hours', 'minutes', 'seconds',
                   'atm pressure', 'atm temperature', 'atm humidity',
                   'T108AP', 'T200Gb', 'T201Gb', 'T202Gb',	'T203Gb',
                   'T300Wa', 'T304Wa', 'T400Wa', 'T405Wa', 'T406Wa',
                   'air flowrate', 'gas flowrate',
                   'F303Wa', 'F401Wa', 'P107Ap',
                   'CO', 'CO2', 'NO', 'O2', 'NO']

fuel = {'CH4': .827, 'C2H6': .039, 'C3H8': .012,
        'CO2': .014, 'N2': .108}

def mixing(fuel, fct):
	result = 0
	for (fl,fraction) in fuel.items():
		result += fraction * fl.fct()
