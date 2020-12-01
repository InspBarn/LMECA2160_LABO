# *-* coding: utf-8 *-*

# Mathematic Developments
import math as mt
import numpy as np

# Figures & plots
from matplotlib import pyplot as plt
try:
    import seaborn as sns
    use_seaborn = True
    ss.set()
except:
    use_seaborn = False

from classCombustible import Combustible
from classMixedFuelTheoretical import MixedFuelTheoretical
from classMixedFuelExperiment import MixedFuelExperiment

from LabResults_donatien import gas_flow_rates, air_flow_rates, T_room, gas_comp as real_gc
gas_flow_rate = gas_flow_rates[0]

disp = True
nfig = 1

#gas_flow_rate = 1.
Ts = 298.15
Ps = 1.01325e5

# Power of the gas supply [kW]
power = 17.5e3
fuel  = {'CH4': .827, 'C2H6': .039, 'C3H8': .012, 'CO2': .014, 'N2': .108}

# Fuel composition as detailed in the .pdf file
mxf = MixedFuelTheoretical(fuel=fuel, Ts=298.15, Ps=1.01325e5) #, AFV_ratio=20./gas_flow_rate)

print("Volumetric Gas flowrate assumption: %.2f [m³/h]" %gas_flow_rate)
print("The volumetric air flowrate for a stoechiometric combustion: %.2f [m³/h]" %(mxf.AFV_stoech*gas_flow_rate))

print("=========== 1 ===========") 
print("  1. The Lower Heating Value (LHV) of the fuel: %.2f [MJ/m³]\n" % (mxf.LHV*mxf.rhom*1e-6))
print("     The Lower Heating Value (LHV) of the fuel: %.2f [MJ/kmol]\n" % (mxf.LHV*1e3))
print("     The theoretical gas flowrate for specified power: %.2f [mol/h]\n" % (3600*power/mxf.LHV))
print("     The theoretical gas flowrate for specified power: %.2f [m³/h]\n" % (3600*power/mxf.LHV/mxf.rhom))
# print("The Lower Heating Value (LHV) of the fuel: %.2f [MJ/kg]\n" % (mxf.LHVM*1E-6))
print("  2. The air-to-fuel ratio at stoechiometry is AFst = %.2f [kg(air)/kg(fuel)]" % mxf.AF_stoech)
print("     The product-to-fuel ratio at stoechiometry is PFst = %.2f [kg(gas)/kg(fuel)]\n" % mxf.PF_stoech)

print("=========== 2 ===========") 
print(mxf.pdt_id)
ost_points = []
for (i,afr) in enumerate(air_flow_rates):
    # Re-set fuel compositon to re-create new mxf with new afr
    # fuel  = {'CH4': .827, 'C2H6': .039, 'C3H8': .012, 'CO2': .014, 'N2': .108}
    # mxf = MixedFuels(fuel=fuel, AFV_ratio=afr/gas_flow_rate)
    AFV_ratio = afr/gas_flow_rates[i]

    mxf._update_AFV_ratio(AFV_ratio,ε=5e-2)
    mxfe = MixedFuelExperiment(real_gc[i], fuel,Ts,Ps,AFV_ratio=AFV_ratio)

    print("  Air flowrate = %.2f" % afr)
    print("    1. The excess-air coefficient λ: %.2f" %mxf.λ)
    print("    2. The flue gas compsition (dry basis):")
    print("       \t\ttheorical \treal")
    ost_points.append({'afr':afr})
    for (comp, val) in mxfe.flue_gas_comp(dry=True).items():
        if comp in ['CO', 'NO']:
            print("       %s:\t%d ppm   \t%d ppm" %(comp, val*1e6, real_gc[i][comp]))
        elif comp in ['CO2', 'O2']:
            print("       %s:\t%.2f %%\t\t%.2f %%" %(comp, val*100, real_gc[i][comp]))
            ost_points[i][comp] = val
        else:
            print("       %s:\t%.2f %%" %(comp, val*100))
    print("    3. Adiabatic combustion temperature: T_ad = %.2f" %mxf.T_ad(T_in=T_room))


if disp:
    fig = plt.figure()
    #nfig += 1
    
    (w,x,y,z) = mxf.get_wxyz()
    P = (0., z/(4.76*z+3.76*(y-2.*x)/4.))
    Q = (1./4.76, 0.)
    S = ((z+y/4.)/(3*(z+y/4.)+7.52*(z+(y-2.*x)/4.)), 0.)
    
    X = [P[0],Q[0],S[0],P[0]]
    Y = [P[1],Q[1],S[1],P[1]]
    
    plt.plot(X,Y,'-o',c='k',lw=1.5)
    plt.xlabel(r"$\left[O_2\right]'$")
    plt.ylabel(r"$\left[CO_2\right]'$")
    plt.title("Ostwald triangle: UCLouvain NG")
    for op in ost_points:
        plt.plot(op['O2'], op['CO2'],'o', label="Air flowrate = %.2f" %op['afr'])
    plt.legend()
    plt.show()

