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

from classCombustible import Combustible
from classMixedFuels import MixedFuels


nfig = 1

# Power of the gas supply [kW]
power = 17.5e3
fuel  = {'CH4': .827, 'C2H6': .039, 'C3H8': .012, 'CO2': .014, 'N2': .108}

# Fuel composition as detailed in the .pdf file
mxf = MixedFuels(fuel=fuel, AFV_ratio=28./20.)

print("Volumetric Gas flowrate assumption : %.2f [m³/h]" %20.)
print("The volumetric air flowrate for a stoechiometric combustion : %.2f [m³/h]" %(mxf.AFV_stoech()*20.))
print("The Lower Heating Value (LHV) of the fuel : %.2f [MJ/m³]\n" % (mxf.LHV_m()*mxf.rho()*1E-6))
#print("The Lower Heating Value (LHV) of the fuel : %.2f [MJ/kg]\n" % (mxf.LHV_m()*1E-6))

print("The air-to-fuel ratio at stoechiometry is AFst = %.2f [kg(air)/kg(fuel)]" % mxf.AF_stoech())
print("The product-to-fuel ratio at stoechiometry is PFst = %.2f [kg(gas)/kg(fuel)]\n" % mxf.PF_stoech())

mxf.excess_air_coefficient()
print("The excess-air coefficient λ : %.2f" %mxf.lambda0)

fig = plt.figure(nfig)
nfig += 1

(w,x,y,z) = mxf.get_xyz()
P = (0., z/(4.76*z+3.76*(y-2.*x)/4.))
Q = (1./4.76, 0.)
S = ((z+y/4.)/(3*(z+y/4.)+7.52*(z+(y-2.*x)/4.)), 0.)

X = [P[0],Q[0],S[0],P[0]]
Y = [P[1],Q[1],S[1],P[1]]

plt.plot(X,Y,'-o',c='k',lw=1.5)
plt.xlabel(r"$\left[O_2\right]'$")
plt.ylabel(r"$\left[CO_2\right]'$")
plt.title("Ostwald's diagram")
plt.show()

