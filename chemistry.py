#
# Cantera .cti input file processor and chemistry solver
#
# The functions and classes in this module process Cantera .cti input
# files and run a simulation. It can be imported as a module, or used
# as a script.
#
# script usage:
# 
# python chemistry.py infile.cti
# 
# Richard West 2010

import sys, os
import math
import pylab, numpy
from scipy.integrate import odeint
import re

# get quantities package from http://pypi.python.org/pypi/quantities
# read about it at http://packages.python.org/quantities/index.html
# I think this may only work with an old version - try "easy_install quantities==0.5b5"
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
    """Was not within valid range for expression"""
    pass
    
class ideal_gas(ctml.ideal_gas):
    """An ideal gas mixture."""
    pass

class state(ctml.state):
    """A state."""
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
        """True if the polynomial is valid at the given temperature, else False"""
        if temperature<self._t[0]: return False
        if temperature>self._t[1]: return False
        return True

    def getThermo(self,T):
        """Get (HeatCapacityOverR, EnthalpyOverRT, EntropyOverR) for a given temperature T
        
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
        """Return the Gibbs Free Energy divided by R (molar gas constant) at temperature T"""
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
    """
    Turn a dictionary  (of concentrations, rates, etc.) into an array AND UNITS.
    
    Gets names (in order) from _speciesnames.
    Returns an array, and a quantities object with the units.
    """
    outArray = pylab.array([inDict[s].simplified for s in _speciesnames])
    # possibly not the fastest way to do it..self.
    units = sum(inDict.values()) / sum(outArray)
    return outArray, units
    
def DictFromArray(inArray, units=None):
    """
    Turn an array (of concentrations, rates, etc.) into a dictionary.
    
    Gets names (in order) from _speciesnames.
    """
    outDict = dict.fromkeys(_speciesnames)
    for i,speciesName in enumerate(_speciesnames):
        value = inArray[i]
        if units:
            value = pq.Quantity(value,units)
        outDict[speciesName] = value
    return outDict

class PropertiesOfSpecies():
    """
    An individual species' properties
    
    >> ps=PropertiesStore()
    >> oxygen = ps['O2(1)'] // oxygen is now a PropertiesOfSpecies instance
    >> oxygen.Radius
    """
    
    # add [property: function] pairs to this dictionary for calculatable fake attributes
    _calculated_properties = dict() 
        
    def __init__(self, properties_dict):
        for key,value in properties_dict.iteritems():
            # try to turn it into a float
            try:
                value = float(value)
            except ValueError:
                pass 
            # insert it into the instance's attribute dictionary
            setattr(self,key,value)
        self.species_name = properties_dict['ChemkinName']
    
    def _getAttributeNames(self):
        """Get a list of extra (fake) attributes. Useful for ipython shell.""" 
        return self.__class__._calculated_properties.keys()
        
    def __getattr__(self,property_name):
        """
        Get a (fake) attribute.
        
        Will only get called if not found in attribute dictionary, or called explicitly.
        Will only call functions stored in the _calculated_properties dictionary of the class. 
        """
        try:
            function = self.__class__._calculated_properties[property_name]
        except KeyError:
            raise AttributeError("Don't know how to calculate %s"%property_name)
        return function(self)
        raise AttributeError # if you haven't already returned
        
    def getMolarVolume(self):
        """Get the molar volume, in m3/mol"""
        molecule_volume =  self.Radius*self.Radius*self.Radius * math.pi * 4/3
        molar_volume = molecule_volume * 6.0221415E23
        return molar_volume / 3.0 #<--- temporary hack!
        # seems to be about 3 times too high.
        # molar volume of undecane = 1/(0.74 g/ml / 156.31 g/mol) = 0.00021122973 m3/mol
        # whereas p['n-C11(2)'].MolarVolume = 0.00063723
    _calculated_properties['MolarVolume']=getMolarVolume  # make a fake attribute
    
    def getPartitionCoefficient(self):
        """
        Get the partition coefficient, K: ratio of solvent to gas concentrations.
        
        K > 1 means concentration is higher in solution than in the gas.
        
        Abraham constants here are for dry decane, taken from from 
        M.H. Abraham et al. / J. Chromatogr. A 1037 (2004) 29-47
        http://dx.doi.org/10.1016/j.chroma.2003.12.004
        """
        # Solvent parameters  here are for dry decane, from from 
        # M.H. Abraham et al. / J. Chromatogr. A 1037 (2004) 29-47
        # http://dx.doi.org/10.1016/j.chroma.2003.12.004
        # if you change them, remember to update the docstring above!
        c = 0.156
        s = 0
        b = 0
        e = -0.143
        l = 0.989
        a = 0
        # Solute parameters are properties of the species
        S=self.AbrahamS
        B=self.AbrahamB
        E=self.AbrahamE
        L=self.AbrahamL
        A=self.AbrahamA
        # Abraham model:
        logK=c + s*S + b*B + e*E + l*L + a*A 
        partition_coefficient = 10**logK
        return partition_coefficient
    # make a fake attribute
    _calculated_properties['PartitionCoefficient']=getPartitionCoefficient 
    
    
    def getDiffusivityInAir(self,Temperature,pressure_in_bar):
        """
        Diffusivity of non-ring HCO compound in air
        
        Use Fuiller, Schettier, and Giddings correlation
        as in Eq. 11-4.1 in Reid, Prausnitz and Sherwood
        
        input: T in K; p in bar;
        output diffusivity d in m2/s
        """        
        wC=12.011; wH=1.008; wO=16;   wair=28.97;
        vC=16.5;   vH=1.98;  vO=5.48; vair=20.1;
        
        vHCO=self.nC*vC+self.nH*vH+self.nO*vO;
        wHCO=self.nC*wC+self.nH*wH+self.nO*wO;
        
        pinv=1./pressure_in_bar;
        
        d=1e-3 * Temperature**1.75 * math.sqrt((wHCO+wair)/(wHCO*wair))*pinv/ \
                (vHCO**.33333+vair**.33333)**2;
        return d*1e-4; #result in m2/s
        

    def getNumberOfCarbons(self):
        """Get the number of Carbon atoms in the molecule"""
        m = re.search('C(\d*)',self.ChemicalFormula)
        if not m: return 0
        if not m.group(1): return 1
        return int(m.group(1))
    _calculated_properties['nC']=getNumberOfCarbons
    def getNumberOfHydrogens(self):
        """Get the number of Carbon atoms in the molecule"""
        m = re.search('H(\d*)',self.ChemicalFormula)
        if not m: return 0
        if not m.group(1): return 1
        return int(m.group(1))
    _calculated_properties['nH']=getNumberOfHydrogens
    def getNumberOfOxygens(self):
        """Get the number of Carbon atoms in the molecule"""
        m = re.search('O(\d*)',self.ChemicalFormula)
        if not m: return 0
        if not m.group(1): return 1
        return int(m.group(1))
    _calculated_properties['nO']=getNumberOfOxygens
    
            
class PropertiesStore():
    """
    A class to store and evaluate Properties of all the Species in the system.
    """
    def __init__(self, resultsDir='RMG_results'):
        self._specs_props = dict()
        self.loadPropertiesFromFile(resultsDir)
        
    def __getattr__(self, property_name):
        """
        Get an array of the property. Length is Nspecies; order is same as chemistry model.
        
        >> s = SpeciesProperties()
        >> s.Radius
        """
        try:
            return self.getPropertyArray(property_name)
        except KeyError:
            raise AttributeError
        
    def __getitem__(self,species_name):
        """Get properties of an individual species."""
        spec_prop = self._specs_props[species_name]
        return spec_prop
        
    def loadPropertiesFromFile(self, resultsDir):
        """
        Load the species properties from the RMG_Solvation_Properties.txt file
        
        This will fill the self.properties dictionary, 
        with dictionaries of properties for each species.
        """
        import csv
        print "Loading species properties from",resultsDir
        propsfilename = os.path.join(resultsDir,'RMG_Solvation_Properties.txt')
        propsfile = open(propsfilename)
        reader = csv.DictReader(propsfile, dialect=csv.excel_tab)
        for spec_prop in reader:
            self._specs_props[spec_prop['ChemkinName']] = PropertiesOfSpecies(spec_prop)
        propsfile.close()
        
    def getSpeciesProperty(self,species_name,property_name):
        """Get the value of a property of a species. General method."""
        if species_name in ['Ar','N2']:
            print "WARNING: Using 'O2(1)' properties for %s because I don't have %s values"%(species_name,species_name)
            species_name='O2(1)'
        try:
            spec_prop = self._specs_props[species_name]
        except KeyError:
            print "Don't have any properties for the species named '%s'."%species_name
            print "These are the species I have:",self._specs_props.keys()
            raise 
        try:
            value = getattr(spec_prop,property_name)
        except AttributeError:
            # print "Don't have the property '%s' for species '%s'."%(property_name,species_name)
            raise
        return value
    
    def getPropertyArray(self,property_name):
        """Get an array of the property. Length is Nspecies; order is same as chemistry model.
        
        Probably quite slow, so wise to store the result."""
        values = numpy.zeros(len(_species))
        for species_index, species_name in enumerate(_speciesnames):
            values[species_index] = self.getSpeciesProperty(species_name,property_name)
        return values

class ChemistrySolver():
    """A chemistry solver."""
    
    def __init__(self, resultsDir='RMG_results'):
        self.loadChemistryModel(resultsDir)
        #: Number of species
        self.Nspecies = len(_species)
        print "Liquid phase has %d species"%self.Nspecies, _speciesnames
        
        self.calculateStoichiometries()
        self.T = 0
        
        self.concentrations = numpy.array([])
        self.properties = PropertiesStore(resultsDir)
    
    def loadChemistryModel(self, resultsDir):
        print "Loading chemistry from",resultsDir
        ctifile = os.path.join(resultsDir,'chemkin','chem.cti')
        base = os.path.basename(ctifile)
        root, ext = os.path.splitext(base)
        ctml.dataset(root)
        if not _species:
            execfile(ctifile)
        else:
            print "Already had chemistry loaded! If you want different chemistry please restart your python shell"
    
    @property
    def speciesnames(self):
        return _speciesnames
    
    def setConcentrations(self, concentrations, zero_others=False):
        """
        Set the concentrations (in mol/m3). Accepts either an array or dictionary.
        
        If passing a dictionary, the keys must be species names, the values their concentrations.
        If zero_others==True then species not in the dictionary will have concentrations set to 0,
        otherwise they will be left alone.
        """
        if type(concentrations)==numpy.ndarray:
            assert concentrations.size==self.Nspecies, "Concentrations array should be Nspecies=%d elements long"%self.Nspecies
            self.concentrations = concentrations
        elif type(concentrations)==dict:
            if zero_others: self.concentrations = numpy.zeros(self.Nspecies)
            assert self.concentrations.size == self.Nspecies, "Concentrations array hasn't been initialised." # try calling with zero_others=True
            for key,value in concentrations.iteritems():
                i = _speciesnames.index(key)
                self.concentrations[i] = value
        else: 
            raise Exception("Expected ethier an array or a dictionary, not a %s"%type(concentrations).__name__)
        print "concentrations", self.concentrations
        
    def getConcentrations(self):
        return self.concentrations
    
    def calculateStoichiometries(self):
        (stoich_reactants,stoich_products,stoich_net) = getStoichiometryArrays()
        self.stoich_reactants = stoich_reactants
        self.stoich_products = stoich_products
        self.stoich_net = stoich_net
        print "stoich_reactants", stoich_reactants
        print "stoich_products", stoich_products
        print "stoich_net", stoich_net
        
    def setTemperature(self, T):
        """
        Set the temperature of the solver, in Kelvin.
        
        If T has changed, recalculates the Rate Coefficients
        """
        
        if T != self.T:
            self.calculateRateCoefficients(T)
            self.T = T
            print "Temperature is now %f K"%T
        else:
            print "Temperature is still %f K"%T
    
    def calculateRateCoefficients(self, T):
        """
        Update the stored forward and reverse rate coefficients.
        
        Call this whenever T changes.
        """
        self.forward_rate_coefficients  = getForwardRateCoefficientsVector(T)
        self.reverse_rate_coefficients  = getReverseRateCoefficientsVector(T)
        print "forward_rate_coefficients", self.forward_rate_coefficients
        print "reverse_rate_coefficients", self.reverse_rate_coefficients
        
    def getRightSideOfODE(self,concentrations,timesteps):
        """Get the function which will be the right side of the ODE"""
        # This is probably the function to speed up, when the time comes for optimization:
        """Get the net rate of creation of all species at a concentration and T."""
        forward_rates = self.forward_rate_coefficients*(concentrations**self.stoich_reactants).prod(1)
        reverse_rates = self.reverse_rate_coefficients*(concentrations**self.stoich_products).prod(1)
        net_rates = forward_rates - reverse_rates
        net_rates_of_creation = numpy.dot(net_rates.T, self.stoich_net)
        return net_rates_of_creation
    
    def getNetRatesOfCreation(self,concentrations):
        """abstracted from getRightSideOfODE
        Get the net rate of creation of all species at a concentration and T."""
        forward_rates = self.forward_rate_coefficients*(concentrations**self.stoich_reactants).prod(1)
        reverse_rates = self.reverse_rate_coefficients*(concentrations**self.stoich_products).prod(1)
        net_rates = forward_rates - reverse_rates
        net_rates_of_creation = numpy.dot(net_rates.T, self.stoich_net)
        return net_rates_of_creation
        
    def solveConcentrationsAfterTime(self, starting_concentrations, reaction_time, temperature=None ):
        """Solve the simulation for a given time starting from a given concentration.
        
        Set the global temperature, if provided (else leave it alone)
        Set the concentrations to starting_concentration.
        Set the time to 0.
        Run the simulation until reaction_time.
        
        Returns the concentrations at the end time."""
        
        if temperature:
            self.setTemperature(T)
        self.setConcentrations(starting_concentrations)
        timesteps = (0, reaction_time)
        concentration_history_array = odeint(self.getRightSideOfODE,starting_concentrations,timesteps)
        return concentration_history_array[-1] # the last row is the final timepoint
        
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
        self.liquidMolarDensity=pq.Quantity(liquidMolarDensity,'mol/m**3') # mol/m3
        self.initialConcentration=self.liquidMolarDensity*initialVolFraction

if __name__ == "__main__":
    import sys, os
    # use different chemistry mechanism if specified on the command line
    if len(sys.argv)>1:
        solver = ChemistrySolver(resultsDir = sys.argv[1])
    else:
        solver = ChemistrySolver()

    # calculate the initial concentrations
    fuel=[
        FuelComponent('n-C11(2)',  0.05,dict(C=11,H=24,O=0),dict(A=6.9722, B=1569.57, C=187.7  ),4945.0),
        FuelComponent('n-C13(3)',  0.19,dict(C=13,H=28,O=0),dict(A=7.00756,B=1690.67, C=174.22 ),4182.0),
        FuelComponent('Mnphtln(4)',0.11,dict(C=11,H=10,O=0),dict(A=7.03592,B=1826.948,C=195.002),3.7e3 ),
        FuelComponent('n-C16(5)',  0.25,dict(C=16,H=34,O=0),dict(A=7.02867,B=1830.51, C=154.45 ),3415.0),
        FuelComponent('C10bnzn(6)',0.12,dict(C=16,H=26,O=0),dict(A=7.8148, B=2396.8,  C=199.5736),2.6e3),
        FuelComponent('n-C19(7)',  0.18,dict(C=19,H=40,O=0),dict(A=7.0153, B=1932.8,  C=137.6  ),2889.0),
        FuelComponent('n-C21(8)',  0.10,dict(C=21,H=44,O=0),dict(A=7.0842, B=2054,    C=120.1  ),2729.0)
    ]
    concs_dict=dict.fromkeys( solver.speciesnames )
    for speciesName in solver.speciesnames:
        concs_dict[speciesName]=pq.Quantity(0.0,'mol/m**3')
    for component in fuel:
        concs_dict[component.name]=component.initialConcentration # mol/m3
    concs_dict['O2(1)']=pq.Quantity(10,'mol/m**3') # haven't a clue
    concentrations,units = ArrayFromDict(concs_dict)
    print "Initial concentrations:", concentrations, units
    
    # Set up solver
    T=430 # kelvin
    solver.setTemperature(T)
    solver.setConcentrations(concentrations)
    
    # set up timesteps
    start=0
    stop=1
    steps=101
    timesteps=pylab.linspace(start,stop,steps)
    
    # solve it here using odeint
    print "Starting to solve it in one go"
    concentration_history_array = odeint(solver.getRightSideOfODE,concentrations,timesteps)
    print "Solved"
    mass_concentrations = concentration_history_array[-1] * solver.properties.MolecularWeight
   # print mass_concentrations
    print 'netRatesofCreation is ', solver.getNetRatesOfCreation(concentrations)
    
    # # plot the graph
    # pylab.semilogy(timesteps,concentration_history_array)
    # pylab.legend(_speciesnames)
    # pylab.show()
    
    # # solve it in the solver, to show the API
    # print "Starting to solve it step by step (in 10 times fewer steps)"
    # time_now = timesteps[0]
    # concentrations_now = concentrations
    # for step in xrange(1,len(timesteps),10):
    #     time = timesteps[step]
    #     concentrations_now = solver.solveConcentrationsAfterTime(concentrations_now, time-time_now )
    #     time_now = time
    #     pylab.semilogy((time_now),concentrations_now.reshape((1,24)), '.')
    # print "Solved"
    
    # mass_concentrations = concentrations * solver.properties.MolecularWeight
    # mass_fractions = mass_concentrations / mass_concentrations.sum()
    
    # gas_phase_concentrations = concentrations / solver.properties.PartitionCoefficient

    
