# *-* coding: utf-8 *-*

import math as mt
import numpy as np

from thermo import Chemical

H2  = Chemical('H2')
O2  = Chemical('O2')
N2  = Chemical('N2')
CO2 = Chemical('CO2')
H2O = Chemical('H2O')

"""
Combustible est un objet représentant un combustible ainsi que sa combustion. Un combustible
est defini par la formule :
                                          Cz Hy Ox
Les paramètres sont :
    - w : indice d'azote
    - x : indice d'oxygene
    - y : indice d'hydrogene
    - z : indice de carbone

Le méthane (CH4) est le combustible defini par defaut.
"""
class Combustible:
    def __init__(self,fuel='CH4', Ts=298.15, Ps=1.0325e5):
        self.chem = Chemical(fuel,Ts,Ps)
        dctn = self.chem.atoms
        self.w = dctn['N'] if 'N' in dctn else 0
        self.x = dctn['O'] if 'O' in dctn else 0
        self.y = dctn['H'] if 'H' in dctn else 0
        self.z = dctn['C'] if 'C' in dctn else 0

    # Return the molar mass in [g/mol]
    @property
    def MW(self):
        return self.chem.MW

    # Return the molar mass in [g/mol]
    @property
    def formula(self):
        return self.chem.formula

    @property
    def atoms(self):
        return self.chem.atoms

    # We assume we have a gas, here (clearly do not want to do more researches !)
    # Return the mass density in [kg/m^3]
    @property
    def rho(self):
        return self.chem.rhog
        # return self.molar_mass / (8.3144621*298.15/1.01325e5)

    # Return the molar density [mol/m^3]
    @property
    def rhom(self):
        return self.chem.rhogm

    # Return the Lower Heating Value (LHV) of the fuel [J/mol]
    @property
    def LHV(self):
        if (self.x == self.y == self.z == 0) \
            or (self.z==1 and self.y==0 and self.x==2) \
            or (self.z==0 and self.y==2 and self.x==1):
            # No Fuel or CO2 or H2O
            return 0
        elif (self.y==2) and (self.x==self.z==0):
        	# Gaseous Fuel // H2
            return 241800
        elif (self.z==self.x==1) and (self.y==0):
        	# Gaseous Fuel // Carbon monoxyde
            return 282400
        else:
            # Gaseous Fuel // Carbon oxydes CnHm
            return self.chem.Hf - self.z*CO2.Hf - self.y/2*H2O.Hf

    # Return the air to fuel ratio of the combustible at stoechiometry
    @property
    def AF_stoech(self):
        if (self.x==self.y==self.z==0):
            return 0
        return (O2.MW + 3.76*N2.MW) * (self.z+(self.y-2*self.x)/4) / self.MW

    # Return the volumetric air to fuel ratio of the combustible at stoechiometry
    @property
    def AFV_stoech(self):
        if (self.x==self.y==self.z==0):
            return 0
        return 4.76 * (self.z + (self.y-2*self.x)/4)

    # Return the product to fuel ratio of the combustible at stoechimotry
    @property
    def PF_stoech(self):
        if (self.x==self.y==self.z==0):
            return 1
        return (self.z*CO2.MW + self.y/2*H2O.MW + 3.76*N2.MW*(self.z+(self.y-2*self.x)/4)) / self.MW
    	#return self.AF_stoech() + 1

    # Return the volumetric product to fuel ratio of the combustible at stoechiometry
    @property
    def PFV_stoech(self):
        return self.z + self.y/2 + 3.76*(self.z + (self.y-2*self.x)/4)

"""
    def balance_equation_to_fuel(eps: float=1E-3):
        left_hand_side = 4.76*self.O2 + (4.76 + 3.76*(y-2*x)/(4*z))*self.CO2 + (2.88 + 3.76*(y-2*x)/(4*z) - 0.88*y/(4*z))*self.CO
        return (left_hand_side - 1 < eps)
"""
