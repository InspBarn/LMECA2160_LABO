# *-* coding: utf-8 *-*

# Read Files
import csv
import pandas as pd

# Mathematic developments
import math as mt
import numpy as np

# Figure and Graphical representation
from matplotlib import pyplot as plt

"""
We first read the file in which results have been saved and save it into a pandas
DataFrame structure.
"""
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
                   'T108AP', 'T200Gb', 'T201Gb', 'T202Gb', 'T203Gb',
                   'T300Wa', 'T304Wa', 'T400Wa', 'T405Wa', 'T406Wa',
                   'air flowrate', 'gas flowrate',
                   'F303Wa', 'F401Wa', 'P107Ap',
                   'CO', 'CO2', 'NO', 'O2', '??']

# Gas Flowrate calibration
results = results[(results['gas flowrate']>0.) & (results['gas flowrate']<5.)]
results.index = range(len(results))
gas_flow_rate = np.mean(results['gas flowrate'])

# Time reset
def reset_time(df):
	df_copy = df.copy()
	time_0 = (df['hours'][0], df['minutes'][0], df['seconds'][0])
	df_copy['hours']   -= time_0[0]
	df_copy['minutes'] -= time_0[1]
	df_copy['seconds'] -= time_0[2]
	return df_copy

def get_flue_gasses(df,fg):
	res = []
	vec = df[fg]
	for el in vec:
		if (not el in res) and (el!=0):
			res.append(el)
	return res

"""
CO  = get_flue_gasses(results,'CO')
CO2 = get_flue_gasses(results,'CO2')
NO  = get_flue_gasses(results,'NO')
O2  = get_flue_gasses(results,'O2')

print((CO,CO2,NO,O2))
"""

"""
Decomposition of the three experimentations as detailed in the .pdf file.
	Expe. 1 -- air flowrate of 20 m^{3} s^{-1}
	Expe. 2 -- air flowrate of 24 m^{3} s^{-1}
	Expe. 3 -- air flowrate of 28 m^{3} s^{-1}
We select a range of 0.4 m^{3} s^{-1} around the air flow rate to allow a certain error.
We want the datas being stabilised to have 
and take the mean of the resulting vector for
	- the air flowrate
	- the gas flowrate
	- T200Gb, T201Gb, T202Gb, T203Gb
	- CO, CO2, NO and O2
We are not interesting in the other datas for our combustion study.
"""
def _select_temperature(df):
	df_copy = df.copy()
	for TGb in ['T200Gb','T201Gb','T202Gb','T203Gb']:
		Tmean = np.mean(df[TGb])
		df_copy = df[(df[TGb] > Tmean*0.95)
		             & (df[TGb] < Tmean*1.05)]
	return df_copy

"""
exp_20 = results[(results['air flowrate'] > 19.8) & (results['air flowrate'] < 20.2)
                      & (results['CO'] != 0) & (results['CO2'] != 0) & (results['O2'] != 0)
                      & (results['O2'] != 0)]
exp_20 = select(exp_20)
"""
"""
experience_24 = results[(results['air flowrate'] > 19.8) & (results['air flowrate'] < 20.2)
                      & (results['CO'] != experience_20['CO'][-1])
                      & (results['CO2'] != experience_20['CO2'][-1])
                      & (results['NO'] != experience_20['NO'][-1])
                      & (results['O2'] != experience_20['O2'][-1])]
experience_24 = select(experience_24)
"""

if __name__=="__main__":
	time = results['hours']*3600 + results['minutes']*60 + results['seconds']
	#time = time - time[0]

	plt.plot(time,results['air flowrate'])

	plt.figure()
	plt.plot(time,results['T201Gb'])
	plt.plot(time,results['T202Gb'])
	plt.plot(time,results['T203Gb'])

	plt.figure()
	plt.plot(time,results['T203Gb'])
	#plt.plot(range(len(omg)),omg['gas flowrate'])
	plt.show()
