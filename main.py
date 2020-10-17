# *-* coding: utf-8 *-*

# Mathematic Developments
import math as mt

# Figures & plots
from matplotlib import pyplot as plt
try:
	import seaborn as sns
	use_seaborn = True
	ss.set()
except:
	use_seaborn = False

import classCombustible as cbt
import classMixedFuels as mxf


# Power of the gas supply [kW]
power = 17.5e3

# Fuel composition as detailed in the .pdf file
fuel = mxf.MixedFuels(fuel = {'CH4': .827, 'C2H6': .039, 'C3H8': .012, 'CO2': .014, 'N2': .108})

print("Volumetric Gas flowrate assumption : 1 [m³/h]")
print("The volumetric air flowrate for a stoechiometric combustion : %.2f [m³/h]" % fuel.AFV_stoech())
print("The Lower Heating Value (LHV) of the fuel : %.2f [MJ/m³]" % (fuel.LHV()*fuel.rho()*1E-6))


