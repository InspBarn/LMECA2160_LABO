# *-* coding: utf-8 *-*

import math as mt
import numpy as np

from scipy.optimize import fsolve

from thermo import Chemical
from classMixedFuel import MixedFuel
from classCombustible import Combustible

H2  = Chemical('H2')
O2  = Chemical('O2')
N2  = Chemical('N2')
CO2 = Chemical('CO2')
H2O = Chemical('H2O')

rho_a = (.21*O2.rho + .79*N2.rho) #/ (8.3144621*298.15/1.01325e5) *1e-3

class MixedFuelTheoretical(MixedFuel):
    def __init__(self, fuel: dict={"CH4": 1.0}, Ts=298.15, Ps=1.01325e5, **kwargs):
        super().__init__(fuel, Ts, Ps)

        if "AF_ratio" in kwargs:
            self._update_AF_ratio(kwargs["AF_ratio"])
        elif "AFV_ratio" in kwargs:
            self._update_AFV_ratio(kwargs["AFV_ratio"])
        else:
            self._update_products()

    #
    def _update_AF_ratio(AF_ratio, ε=0):
        self.AF_ratio = AF_ratio
        self._update_ractants()
        self._update_products(ε)

    # 
    def _update_AFV_ratio(self, AFV_ratio, ε=0):
        self.AF_ratio = AFV_ratio * rho_a/self.rho
        self._update_reactants()
        self._update_products(ε)

    # Update the flue gases composition, with coefficients in stoechiometry with x, y, z
    def _update_products(self, ε=0):
        w,x,y,z = self.get_wxyz()
        if self.λ*(1+ε) < 1:
            # Case 1 : Rich
            k = 2*(1-self.λ)*(4*z+y-2*x)/(4*z+y)
        elif 1 <= self.λ*(1-ε):
            # Case 2 : Poor
            k = 0
        else:
            # Case 3 : Mixed
            k = ((1-self.λ*(1-ε))**2)/(2*ε*self.λ) * (4*z+y-2*x)/(4*z+y)

        self.pdt_frac[self.pdt_id['CO2']] = (1-k)*z
        self.pdt_frac[self.pdt_id['CO']]  = k*z
        self.pdt_frac[self.pdt_id['H2']]  = k*y/4
        self.pdt_frac[self.pdt_id['H2O']] = (1-k/2)*y/2
        self.pdt_frac[self.pdt_id['O2']]  = (self.λ-1)*(z+(y-2*x)/4) + k/2*(z+y/4)
        self.pdt_frac[self.pdt_id['N2']]  = 3.76*self.λ*(z+(y-2*x)/4) + w/2
