# *-* coding: utf-8 *-*

import math as mt
import numpy as np

from scipy.optimize import fsolve

from thermo import Chemical
from classCombustible import Combustible

H2  = Chemical('H2')
O2  = Chemical('O2')
N2  = Chemical('N2')
CO2 = Chemical('CO2')
H2O = Chemical('H2O')

rho_a = (.21*O2.rho + .79*N2.rho) #/ (8.3144621*298.15/1.01325e5) *1e-3

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
    def __init__(self, fuel: dict={"CH4": 1.0}, Ts=298.15, Ps=1.01325e5, **kwargs):
        self.Ts = Ts
        self.Ps = Ps

        # Fuel Decomposition
        self.fuel_id = {fid:i for (i,fid) in enumerate(fuel.keys())}
        self.fuel_n  = len(self.fuel_id)
        self.fuel_type = []
        self.fuel_frac = []
        for (ft,ff) in fuel.items():
            self.fuel_type.append(Combustible(ft,Ts,Ps))
            self.fuel_frac.append(ff)
        self.fuel_frac = np.array(self.fuel_frac)
        self.fuel_frac /= self.fuel_frac.sum()

        ##
        if 'AF_ratio' in kwargs:
            self.AF_ratio = kwargs["AF_ratio"]
        elif 'AFV_ratio' in kwargs:
            self.AF_ratio = kwargs["AFV_ratio"]  * rho_a/self.rho
        else:
            self.AF_ratio = self.AF_stoech

        ## Reactants Decomposition
        self.rct_id = {ft:i for (i,ft) in enumerate(fuel)}
        self.rct_n  = len(fuel)
        if not('N2' in fuel):
            self.rct_id['N2'] = self.rct_n
            self.rct_n += 1
        if not('O2' in fuel):
            self.rct_id['O2'] = self.rct_n
            self.rct_n += 1

        w,x,y,z = self.get_wxyz()
        self.rct_chem = [Chemical(rid,Ts,Ps) for rid in self.rct_id.keys()]
        self.rct_frac = np.zeros(self.rct_n)
        for rid,i in self.rct_id.items():
            if rid=='O2':
                self.rct_frac[i] = self.λ*(z+(y-2*x)/4)
                if 'O2' in fuel: self.rct_frac[i] += fuel['O2']
            elif rid=='N2':
                self.rct_frac[i] = 3.76*self.λ*(z+(y-2*x)/4)
                if 'N2' in fuel: self.rct_frac[i] += fuel['N2']
            else:
                self.rct_frac[i] = fuel[rid]

        ## Products Creation
        self.pdt_id = {'CO2':0, 'CO':1, 'H2':2, 'H2O':3, 'O2':4, 'N2':5}
        self.pdt_n  = len(self.pdt_id)
        self.pdt_chem = [Chemical(pid,Ts,Ps) for pid in self.pdt_id.keys()]
        self.pdt_frac = np.zeros(self.pdt_n)
        self._update_products()

        """
        glob = sum(self.reactant.values())
        for (chem,comp) in self.reactant.items():
            self.reactant[chem] = comp / glob

        glob = sum(self.product.values())
        for (chem,comp) in self.product.items():
            self.product[chem] = comp / glob
        """

    def get_wxyz(self):
        w,x,y,z = 0,0,0,0
        for (ftype, ffraction) in zip(self.fuel_type, self.fuel_frac):
            w += ftype.w * ffraction
            x += ftype.x * ffraction
            y += ftype.y * ffraction
            z += ftype.z * ffraction
        return (w,x,y,z)

    # Return the molar mass in [g/mol]
    @property
    def MW(self):
        result = 0.0
        for (ftype, ffraction) in zip(self.fuel_type, self.fuel_frac):
            result += ftype.MW * ffraction
        return result

    # Return the mass density of the fuel [kg/m^3]
    @property
    def rho(self):
        if not hasattr(self, 'fuel_fraction_m'):
            self._cpt_Mfraction()
        result = 0.0
        for (ftype, ffraction) in zip(self.fuel_type, self.fuel_fraction_m):
            result += (ftype.rho)**(-1) * ffraction
        return (result)**(-1)

    # Return the molar density of the fuel [mol/m^3]
    @property
    def rhom(self):
        result = 0.0
        for (ftype, ffraction) in zip(self.fuel_type, self.fuel_frac):
            result += (ftype.rhom)**(-1) * ffraction
        return (result)**(-1)

    # Compute the mass fraction of fuel components and create a new attribute
    def _cpt_Mfraction(self):
        result = []
        for (ftype, ffraction) in zip(self.fuel_type, self.fuel_frac):
            result.append(ffraction * ftype.MW/self.MW)
        self.fuel_fraction_m = result
        return result

    # Compute the volumetric fraction of fuel components and create a new attribute
    def _cpt_Vfraction(self):
        result = []
        for (ftype, ffraction) in zip(self.fuel_type, self.fuel_frac):
            result.append(ffraction / (ftype.rhom/self.rhom))
        self.fuel_fraction_v = result
        return result

    # Return the fuel LHV [J/mol]
    @property
    def LHV(self):
        result = 0.0
        for (ftype, ffraction) in zip(self.fuel_type, self.fuel_frac):
            result += ftype.LHV * ffraction
        return result

    # Return the fuel LHV [J/kg]
    @property
    def LHVM(self):
        return self.LHV / (self.MW*1e-3)

    # Return the air to fuel ratio at stoechiometry
    @property
    def AF_stoech(self):
        if not hasattr(self, 'fuel_fraction_m'):
            self._cpt_Mfraction()
        result = 0.0
        for (ftype, ffraction) in zip(self.fuel_type, self.fuel_fraction_m):
            result += ftype.AF_stoech * ffraction
        return result

    # Return the volumetric air to fuel ratio at stoechiometry
    @property
    def AFV_stoech(self):
        if not hasattr(self, 'fuel_fraction_v'):
            self._cpt_Vfraction()
        result = 0.0
        for (ftype, ffraction) in zip(self.fuel_type, self.fuel_fraction_v):
            result += ftype.AFV_stoech * ffraction
        return result

    # Return the product to fuel ratio at stoechiometry
    @property
    def PF_stoech(self):
        if not hasattr(self, 'fuel_fraction_m'):
            self._cpt_Mfraction()
        result = 0.0
        for (ftype, ffraction) in zip(self.fuel_type, self.fuel_fraction_m):
            result += ftype.PF_stoech * ffraction
        return result

    # Return the volumetric product to fuel ratio at stoechiometry
    @property
    def PFV_stoech(self):
        if not hasattr(self, 'fuel_fraction_v'):
            self._cpt_Vfraction()
        result = 0.0
        for (ftype, ffraction) in zip(self.fuel_type, self.fuel_fraction_v):
            result += ftype.PFV_stoech * ffraction
        return result

    # Return the excess air coefficient of the combustion
    @property
    def λ(self):
        return self.AF_ratio / self.AF_stoech

    def optimal_excess_air_coefficient(self):
        pass

    """
    # Update the CO coefficient
    @property
    def carbon_monoxide_coefficient(self):
        self.k = self.product['CO'] / (self.product['CO'] + self.product['CO2'])
    """

    # Return the equivalence ratio
    @property
    def equivalence_ratio(self):
        return (self.λ)**(-1)

    # 
    def _update_AFV_ratio(self, AFV_ratio):
        self.AF_ratio = AFV_ratio * rho_a/self.rho
        w,x,y,z = self.get_wxyz()

        #self.rct_frac[self.rct_id['O2']] = self.AF_ratio*O2.MW/(O2.MW+3.76*N2.MW)
        self.rct_frac[self.rct_id['O2']] = self.λ*(z+(y-2*x)/4)
        if 'O2' in self.fuel_id:
            self.rct_frac[self.rct_id['O2']] += self.fuel_frac[self.fuel_id['O2']]

        #self.rct_frac[self.rct_id['N2']] = self.AF_ratio*(3.76*N2.MW)/(O2.MW+3.76*N2.MW)
        self.rct_frac[self.rct_id['N2']] = 3.76*self.λ*(z+(y-2*x)/4)
        if 'N2' in self.fuel_id:
            self.rct_frac[self.rct_id['N2']] += self.fuel_frac[self.fuel_id['N2']]

        self._update_products()


    # Return the flue gases composition, with coefficients in stoechiometry with x, y, z
    def _update_products(self,ε=0):
        w,x,y,z = self.get_wxyz()

        if self.λ*(1+ε) < 1:     # Case 1 : Rich
            k = 2*(1-self.λ)*(4*z+y-2*x)/(4*z+y)
        elif 1 < self.λ*(1-ε):   # Case 2 : Poor
            k = 0
        else:                    # Case 3 : Mixed
            k = ((1-self.λ*(1-ε))**2)/(2*ε*self.λ) * (4*z+y-2*x)/(4*z+y)

        self.pdt_frac[self.pdt_id['CO2']] = (1-k)*z
        self.pdt_frac[self.pdt_id['CO']]  = k*z
        self.pdt_frac[self.pdt_id['H2']]  = k*y/4
        self.pdt_frac[self.pdt_id['H2O']] = (1-k/2)*y/2
        self.pdt_frac[self.pdt_id['O2']]  = (self.λ-1)*(z+(y-2*x)/4) + k/2*(z+y/4)
        self.pdt_frac[self.pdt_id['N2']]  = 3.76*self.λ*(z+(y-2*x)/4) + w/2

    """
    Return the flue gases composition (sum = 1)
    """
    def flue_gas_comp(self,ε=0,dry=False):
        self._update_products(ε)
        total = sum(self.pdt_frac)
        if dry:
            total -= self.pdt_frac[self.pdt_id['H2O']]
        fgc = {}
        for (pid,i) in self.pdt_id.items():
            if not(dry & (pid=='H2O')):
                fgc[pid] = self.pdt_frac[i] / total
        return fgc

    """
    Return the adiabatic temperature of the combustion
    """
    def T_ad(self, T_in=300):
        self._update_products()

        #"""
        HCGr = self.rct_frac @ np.array(
            [chem.HeatCapacityGas.T_dependent_property_integral(self.Ts, T_in)/(T_in-self.Ts)
            for chem in self.rct_chem])

        HCGp = lambda Tad: self.pdt_frac @ np.array(
            [chem.HeatCapacityGas.T_dependent_property_integral(self.Ts, Tad)/(Tad-self.Ts)
            for chem in self.pdt_chem])

        func = lambda Tad: (Tad-self.Ts)*HCGp(Tad) - (self.LHV + (T_in-self.Ts)*HCGr)
        T_ad = fsolve(func, 1000)[0]

        """
        sum_reac = 0
        for (frac,chem) in zip(self.rct_frac,self.rct_chem):
            sum_reac += frac * chem.HeatCapacityGas.T_dependent_property_integral(self.Ts, T_in)/(T_in-self.Ts)

        n,nmax = 0,100
        T_ad = 500
        T_ad_p = 0
        #breakpoint()
        while (mt.fabs(T_ad-T_ad_p)>1e-2) and (n<nmax):
            T_ad_p = T_ad

            sum_prod = 0
            for (frac,chem) in zip(self.rct_frac,self.rct_chem):
                sum_prod += frac * chem.HeatCapacityGas.T_dependent_property_integral(self.Ts, T_ad)/(T_ad-self.Ts)

            T_ad = self.Ts + (self.LHV + (T_in-self.Ts)*sum_reac)/sum_prod
            n += 1

            #print(T_ad) # remove when working

        if n==nmax:
            print("Too much iteration: exiting with T_ad = %.2f" %T_ad)
        #"""

        return T_ad
