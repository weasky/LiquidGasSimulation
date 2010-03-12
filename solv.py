##
# BASED ON: ctml.py
#
# Cantera .cti input file processor
#
# The functions and classes in this module process Cantera .cti input
# files and run a simulation. It can be imported as a module, or used
# as a script.
#
# script usage:
# 
# python solv.py infile.cti
# 
# This will do something
#
# Richard West 2009 

import math
import pylab, numpy
from scipy.integrate import odeint

# get quantities package from http://pypi.python.org/pypi/quantities
# read about it at http://packages.python.org/quantities/index.html
import quantities as pq
pq.UnitQuantity('kilocalories', pq.cal*1e3, symbol='kcal')
pq.UnitQuantity('kilojoules', pq.J*1e3, symbol='kJ')
pq.UnitQuantity('kilomoles', pq.mol*1e3, symbol='kmol')


# import all the classes which write XML
# and all the constant, functions, etc.
#reload(ctml)
#from ctml_writer import *
#import ctml_writer
# then we can modify the ones we want to do more
import ctml_writer as ctml
from ctml_writer import units, OneAtm

_uConc=pq.Quantity(1,ctml._umol) / pq.Quantity(1,ctml._ulen)**3
_uTime=pq.Quantity(1,ctml._utime)

_species=ctml._species
_reactions=ctml._reactions
_speciesnames=ctml._speciesnames


class outsideValidRangeError(Exception):
    """was not within valid range for expression"""
    pass
    
class ideal_gas(ctml.ideal_gas):
    pass

class state(ctml.state):
    pass
    
    
class species(ctml.species):
    """ species derived from ctml_writer.species.
    species(name,atoms,note,thermo,transport,charge)
    """
    def getCorrectNasaPoly(self,T):
        for nasapoly in self._thermo:
            if nasapoly.ValidTemperature(T): return nasapoly
        return None
                    
    def getGibbsFreeEnergy(self,T):
        """getGibbsFreeEnergy(T) returns  FreeEnergyOverR at temperature T"""
        return self.getCorrectNasaPoly(T).getGibbsFreeEnergy(T)
    
    def getThermo(self,T):
        thermo=self.getCorrectNasaPoly(T).getThermo(T)
        return thermo

class NASA(ctml.NASA):
    """NASA polynomial representation of thermo"""
    def ValidTemperature(self,temperature):
        if temperature<self._t[0]: return False
        if temperature>self._t[1]: return False
        return True

    def getThermo(self,T):
        """ getThermo(T) returns 
            (HeatCapacityOverR, EnthalpyOverRT, EntropyOverR) 
        for a given temperature T
        Raises outsideValidRangeError exception if T is not within 
        range of polynomial"""
        if not self.ValidTemperature(T): raise outsideValidRangeError
        if self._pref > 0.0: raise Exception("not sure what to do \
            with customised standard state pressure")
        # http://www.me.berkeley.edu/gri-mech/data/nasa_plnm.html
        # Cp/R = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4
        # H/RT = a1 + a2 T /2 + a3 T^2 /3 + a4 T^3 /4 + a5 T^4 /5 + a6/T
        # S/R  = a1 lnT + a2 T + a3 T^2 /2 + a4 T^3 /3 + a5 T^4 /4 + a7
        c=self._coeffs
        HeatCapacityOverR=0.0
        EnthalpyOverRT=0.0
        EntropyOverR=0.0
        for i in range(5):
            HeatCapacityOverR += c[i] * T**i
            EnthalpyOverRT += c[i] * T**i / (i+1)
        EnthalpyOverRT += c[5] / T
        EntropyOverR = ( c[0]*math.log(T) + c[1]*T + c[2]*T*T/2 +
                            c[3]*T*T*T/3 + c[4]*T*T*T*T/4 + c[6] )
        return(HeatCapacityOverR,EnthalpyOverRT,EntropyOverR)
        
    def getGibbsFreeEnergy(self,T):
        """getGibbsFreeEnergy(T) returns  FreeEnergyOverR at temperature T"""
        (HeatCapacityOverR,EnthalpyOverRT,EntropyOverR) = self.getThermo(T)
        FreeEnergyOverRT = EnthalpyOverRT - EntropyOverR
        return FreeEnergyOverRT*T
        

# Classes should be named with CamelCase but need to stick with Cantera .ctml convention, hence lowercase:
class reaction(ctml.reaction):
    """A chemical reaction."""
    
    def getForwardRateCoefficient(self,T):
        """Get the forward rate coefficient at a given temperature 
        
        assuming kf = [A, n, E] and kf=A (T/1K)^n exp(-E/RT)"""
        
        kf=self._kf
        forwardRateCoefficient=kf[0] / _uConc**(self.getReactantNu()-1) / _uTime
        if kf[1]:
            forwardRateCoefficient *= T**kf[1]
        #assert ctml._ue=='kcal/mol', 'unit conversion not implemented'
        #R=0.0019872065 # kcal/mol/K
        R=pq.constants.R 
        Ea = pq.Quantity(kf[2],ctml._ue)
        Temp = pq.Quantity(T,'K')
        forwardRateCoefficient*=math.exp(-Ea/(R*Temp))
        # TODO: this is where we should correct for diffusion limit (though it may depend on reverse also!)
        print "%s forward k = %s"%(self._e,forwardRateCoefficient)
        return forwardRateCoefficient
        
    def getReverseRateCoefficient(self,T):
        """Get the reverse rate coefficient at a given temperature"""
        if not self.rev: 
            print "%s is irreversible"%(self._e)
            return pq.Quantity(0,"")
        forwardRateCoefficient=self.getForwardRateCoefficient(T)
        reverseRateCoefficient=forwardRateCoefficient/self.getEquilibriumConstant(T)
        print "%s reverse k = %s"%(self._e,reverseRateCoefficient)
        return reverseRateCoefficient
        
    def getDeltaGoverR(self,T):
        """Get the change in Gibbs free energy *over R* at a given T"""
        deltaGOverR=0
        for speciesName,order in self._p.items(): # add products
            deltaGOverR += order* getSpeciesByName(speciesName).getGibbsFreeEnergy(T)
        for speciesName,order in self._r.items(): # subtract reactants
            deltaGOverR -= order* getSpeciesByName(speciesName).getGibbsFreeEnergy(T)       
        return deltaGOverR
        
    def getEquilibriumConstant(self,T):
        """returns the equilibrium constant at a given temperature:  
        Keq = exp(-DeltaG/RT) """
        
        #if not sum(self._p.values())==sum(self._rxnorder.values()):
        #    print "getEquilibriumConstant currently assumes forward and reverse reactions have same reaction order"
        
        dGoverR=self.getDeltaGoverR(T)
        Keq= math.exp(-dGoverR/T) 
        
        Kc=Keq * _uConc ** self.getDeltaNu()
        if self.getDeltaNu(): print "Warning: assuming standard state and real concentrations equal" # or something like that!
        # TODO: the units are right, but it's not _uConc that you should be multiplying by
        
        return Kc
    
    def getReactantNu(self):
        """Get the stoichiometry in the forwards direction.
        
        i.e. the number of reactant molecules.
        uses self._reactantNu to cache the answer"""
        if hasattr(self,'_reactantNu'): return self._reactantNu
        reactantNu=0
        for speciesName,order in self._r.items(): 
            reactantNu += order
        self._reactantNu=reactantNu
        return reactantNu
    def getProductNu(self):
        """Get the stoichiometry in the reverse direction.
        
        i.e. the number of product molecules.
        uses self._productNu to cache the answer"""
        if hasattr(self,'_productNu'): return self._productNu
        productNu=0
        for speciesName,order in self._p.items():
            productNu += order
        self._productNu=productNu
        return productNu
    def getDeltaNu(self):
        """Get the change in stoichiometry of the reaction, delta Nu.
        
        deltaNu= productNu - reactantNu
        uses self._deltaNu to cache the answer"""
        if hasattr(self,'_deltaNu'): return self._deltaNu
        self._deltaNu=self.getProductNu()-self.getReactantNu()
        return self._deltaNu
        
    def getStoichiometryReactantsRow(self):
        """Get the stoichiometry of each species as a reactant"""
        row = numpy.zeros(len(_species))
        for species_name,order in self._r.items(): 
            row[_speciesnames.index(species_name)] = order
        return row
    def getStoichiometryProductsRow(self):
        """Get the stoichiometry of each species as a product"""
        row = numpy.zeros(len(_species))
        for species_name,order in self._p.items(): 
            row[_speciesnames.index(species_name)] = order
        return row
    def getStoichiometryNetRow(self):
        """Get the net stoichiometry of each species (products - reactants)"""
        row = self.getStoichiometryProductsRow() - self.getStoichiometryReactantsRow()
        return row
        
class FuelComponent():
    """Shouldn't really be part of the solv package. 
    specific to the fuel model."""
    def __str__(self):
        return "<Species %s>"%(self.name)
    def __init__(self, name="",
            initialVolFraction=0,
            composition=dict(C=0,H=0,O=0), 
            Antoine=dict(A=0,B=0,C=0),
            liquidMolarDensity=3500 ):
        self.name=name
        self.initialVolFraction=initialVolFraction
        self.composition=composition
        self.Antoine=Antoine
        self.liquidMolarDensity=liquidMolarDensity * pq.mol / pq.m**3  # mol/m3
        self.initialConcentration=self.liquidMolarDensity*initialVolFraction
        
        
## Functions with a global scope:
# TODO: put these in a class?
def getStoichiometryArrays():
    """
    Get arrays of stoichiometric coefficients for reactants, products, and net change.
    
    Returns three arrays: stoich_reactants, stoich_products, stoich_net
    These each have the size (Nreactions x Nspecies)
    """
    stoich_reactants = numpy.zeros((len(_reactions),len(_species)))
    stoich_products = numpy.zeros_like(stoich_reactants)
    stoich_net = numpy.zeros_like(stoich_reactants)
    
    for rn, r in enumerate(_reactions):
        stoich_reactants[rn] = r.getStoichiometryReactantsRow()
        stoich_products[rn] = r.getStoichiometryProductsRow()
        stoich_net[rn] = r.getStoichiometryNetRow()
    return stoich_reactants,stoich_products,stoich_net
    
def getForwardRateCoefficientsVector(T):
    """Get the vector of forward rate coefficients."""
    forward_rate_coefficients = numpy.zeros(len(_reactions))
    for rn, r in enumerate(_reactions):
        forward_rate_coefficients[rn] = r.getForwardRateCoefficient(T).simplified
    return forward_rate_coefficients
        
def getReverseRateCoefficientsVector(T):
    """Get the vector of reverse rate coefficients."""
    reverse_rate_coefficients = numpy.zeros(len(_reactions))
    for rn, r in enumerate(_reactions):
        reverse_rate_coefficients[rn] = r.getReverseRateCoefficient(T).simplified
    return reverse_rate_coefficients
    
def getSpeciesByName(name):
    """Select a species by its name."""
    for s in _species:
        if s._name==name:
            return s
    return None

def ArrayFromDict(inDict):
    """Turn a dictionary  (of concentrations, rates, etc.) into an array AND UNITS.
    
    Gets names (in order) from _speciesnames.
    Returns an array, and a quantities object with the units."""
    outArray = pylab.array([inDict[s].simplified for s in _speciesnames])
    # possibly not the fastest way to do it..self.
    units = sum(inDict.values()) / sum(outArray)
    return outArray, units
    
def DictFromArray(inArray, units=None):
    """Turns an array (of concentrations, rates, etc.) into a dictionary.
    
    Gets names (in order) from _speciesnames"""

    outDict = dict.fromkeys(_speciesnames)
    for i,speciesName in enumerate(_speciesnames):
        value = inArray[i]
        if units:
            value = pq.Quantity(value,units)
        outDict[speciesName] = value
    return outDict


        
if __name__ == "__main__":
    #reload(ctml)
    import sys, os
    if len(sys.argv)>1:
        file = sys.argv[1]
    else:
        print "using default file, chem.cti"
        file = 'chem.cti' # default input file
    base = os.path.basename(file)
    root, ext = os.path.splitext(base)
    ctml.dataset(root)
    if not _species:
        execfile(file)
    
    # set up the concentrations
    fuel=[
        FuelComponent('n-undecane(2)',0.05,dict(C=11,H=24,O=0),dict(A=6.9722,B=1569.57,C=187.7 ),4945.0),
        FuelComponent('n-tridecane(3)',0.19,dict(C=13,H=28,O=0),dict(A=7.00756,B=1690.67,C=174.22),4182.0),
        FuelComponent('SPC(4)',0.11,dict(C=11,H=10,O=0),dict(A=7.03592,B=1826.948,C=195.002),3.7e3),
        FuelComponent('n-hexadecane(5)',0.25,dict(C=16,H=34,O=0),dict(A=7.02867,B=1830.51,C=154.45),3415.0),
        FuelComponent('SPC(6)',0.12,dict(C=16,H=26,O=0),dict(A=7.8148,B=2396.8,C=199.5736),2.6e3),
        FuelComponent('n-nonadecane(7)',0.18,dict(C=19,H=40,O=0),dict(A=7.0153,B=1932.8,C=137.6 ),2889.0),
        FuelComponent('n-heneicosane(8)',0.10,dict(C=21,H=44,O=0),dict(A=7.0842,B=2054,C=120.1),2729.0)
    ]
    concs=dict.fromkeys(_speciesnames)
    for speciesName in _speciesnames:
        concs[speciesName]=pq.Quantity(0.0,'mol/m**3')
    for component in fuel:
        concs[component.name]=component.initialConcentration # mol/m3
    concs['O2(1)']=pq.Quantity(10,'mol/m**3') # haven't a clue
    
    concentrations,units = ArrayFromDict(concs)
    print "concentrations:", concentrations, units
    
    # set up timesteps
    start=0
    stop=100
    steps=1000
    timesteps=pylab.linspace(start,stop,steps)
    
    T=430 # kelvin
    
    (stoich_reactants,stoich_products,stoich_net) = getStoichiometryArrays()
    print "stoich_reactants", stoich_reactants
    print "stoich_products", stoich_products
    print "stoich_net", stoich_net
    
    forward_rate_coefficients  = getForwardRateCoefficientsVector(T)
    print "forward_rate_coefficients", forward_rate_coefficients
    reverse_rate_coefficients  = getReverseRateCoefficientsVector(T)
    print "reverse_rate_coefficients", reverse_rate_coefficients
   # log_concentrations = numpy.log10([c for c in concentrations])
    
    def RightSideOfODE(concentrations, time):
        """Get the net rate of creation of all species at a concentration and T."""
        forward_rates = forward_rate_coefficients*(concentrations**stoich_reactants).prod(1)
        reverse_rates = reverse_rate_coefficients*(concentrations**stoich_products).prod(1)
        net_rates = forward_rates - reverse_rates
        net_rates_of_creation = numpy.dot(net_rates.T, stoich_net)
        return net_rates_of_creation
    
    concentration_history_array = odeint(RightSideOfODE,concentrations,timesteps)
    
    pylab.semilogy(timesteps,concentration_history_array)
    pylab.show()
    pylab.legend(_speciesnames)
    

