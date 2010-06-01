#
# load, store and calculate properties of species in the simulation
# 
# Richard West 2010

import sys, os
import math
import pylab, numpy
import re
import types
    
def getSpeciesByName(name):
    """Select a species by its name."""
    for s in _species:
        if s._name==name:
            return s
    return None

class PropertiesOfSpecies(object):
    """
    An individual species' properties
    
    >> ps=PropertiesStore()
    >> oxygen = ps['O2(1)'] // oxygen is now a PropertiesOfSpecies instance
    >> oxygen.Radius
    """

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

    def getMolarVolume(self):
        """Get the molar volume, in m3/mol"""
        molecule_volume =  self.Radius*self.Radius*self.Radius * math.pi * 4/3
        molar_volume = molecule_volume * 6.0221415E23
        return molar_volume / 3.0 #<--- temporary hack!
        # seems to be about 3 times too high.
        # molar volume of undecane = 1/(0.74 g/ml / 156.31 g/mol) = 0.00021122973 m3/mol
        # whereas p['n-C11(2)'].MolarVolume = 0.00063723
    MolarVolume = property(getMolarVolume)  # make a fake attribute
    
    def getPartitionCoefficient298(self):
        """
        Get the solution/gas partition coefficient, K, at 298K
        
        K = ratio of solvent to gas concentrations.
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
    PartitionCoefficient298 = property(getPartitionCoefficient298)
    
    def PartitionCoefficientT(self, T):
        """
        Get the solution/gas partition coefficient at the given temperature.
        
        K: ratio of solvent to gas concentrations.
        K > 1 means concentration is higher in solution than in the gas.
        See getPartitionCoefficient298 for more details.
        """
        
        deltaH0, deltaS0 = self.SolvationThermochemistry
        deltaG0 = deltaH0 - T * deltaS0
        lnK = deltaG0 / (-8.314 * T )
        partition_coefficient = math.exp(lnK)
        return partition_coefficient
    
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
    nC = property(getNumberOfCarbons)
    def getNumberOfHydrogens(self):
        """Get the number of Carbon atoms in the molecule"""
        m = re.search('H(\d*)',self.ChemicalFormula)
        if not m: return 0
        if not m.group(1): return 1
        return int(m.group(1))
    nH = property(getNumberOfHydrogens)
    def getNumberOfOxygens(self):
        """Get the number of Carbon atoms in the molecule"""
        m = re.search('O(\d*)',self.ChemicalFormula)
        if not m: return 0
        if not m.group(1): return 1
        return int(m.group(1))
    nO = property(getNumberOfOxygens)
    
    def getSolvationThermochemistryMintz(self):
        """
        Get the solvation enthalpy and entropy. Free Energy from Abraham; Enthlapy from Mintz.
        
        In this model neither is temperature-dependent.
        Returns a tuple: (DHsolv, DSsolv)
        
        DHsolv is in J/mol
        DSsolv is in J/mol/K
        """
        deltaG0 = self.getSolvationFreeEnergy298()
        deltaH0 = self.getSolvationEnthalpyMintz()
        T = 298 # standard state temperature
        deltaS0 = (deltaH0 - deltaG0 ) / T
        return deltaH0, deltaS0
        
    def getSolvationThermochemistryPierotti(self):
        """
        Get the solvation enthalpy and entropy. Free Eneregy from Abraham; Entropy from Pierotti
        
        In this model neither is temperature-dependent.
        Returns a tuple: (DHsolv, DSsolv)
        
        DHsolv is in J/mol
        DSsolv is in J/mol/K
        """
        deltaG0 = self.getSolvationFreeEnergy298()
        deltaS0 = self.getSolvationEntropyPierotti()
        T = 298 # standard state temperature
        deltaH0 = deltaG0 + (T*deltaS0)
        return deltaH0, deltaS0

    def getSolvationFreeEnergy298(self):
        """
        Get the change in Gibbs free energy due to solvation at 298K.
        
        Returns Delta G at 298K in J/mol.
        """
        lnK = math.log( self.getPartitionCoefficient298() ) 
        deltaG0 = -8.314 * 298 * lnK;
        return deltaG0
    
    def getSolvationEntropyPierotti(self):
        """
        Get the solvation entropy, from Pierotti's scaled particle theory.
        
        In this model it's not temperature-dependent.
        
        Returns $Delta S_{solv}$ in J/mol/K
        
        Based on code from RMG-Java, based on Rob Ashcraft's thesis around page 60,
        based on THE SOLUBILITY OF GASES IN LIQUIDS Robert A. Pierotti (1963)
        Journal of Physical Chemistry 67 (9) p. 1840-1845
        http://dx.doi.org/10.1021/j100803a024
        """
        r_solute = self.Radius # should be in metres
        r_solvent = 3.498E-10  # Manually assigned solvent radius [=] meter Calculated using Connolly solvent excluded volume from Chem3dPro
        r_cavity = r_solute + r_solvent;   # Cavity radius [=] meter
        rho = 3.09E27   # number density of solvent [=] molecules/Angstrom^3   Value here is for decane using density =0.73 g/cm3
        parameter_y = 4.1887902*rho* r_solvent*r_solvent*r_solvent #  Parameter y from Ashcraft Thesis Refer pg no. 60. (4/3)*pi*rho*r^3
        parameter_ymod = parameter_y/(1-parameter_y) # parameter_ymod= y/(1-y) Defined for convenience
        R=8.314 # Gas constant units J/mol K
        # Definitions of K0, K1 and K2 correspond to those for K0', K1' and K2' respectively from Ashcraft's Thesis (-d/dT of K0,K1,K2)
        K0 = -R*(-math.log(1-parameter_y)+(4.5*parameter_ymod*parameter_ymod))
        K1 = (R*0.5/r_solvent)*((6*parameter_ymod)+(18*parameter_ymod*parameter_ymod))
        K2 = -(R*0.25/(r_solvent*r_solvent))*((12*parameter_ymod)+(18*parameter_ymod*parameter_ymod))
        #Basic definition of entropy change of solvation from Ashcfrat's Thesis
        deltaS0 = K0+(K1*r_cavity)+(K2*r_cavity*r_cavity)
        return deltaS0
    
    def getSolvationEnthalpyPierotti(self):
        """Get the solvation enthlapy, in J/mol, from Abraham's G and Pierotti's S"""
        deltaH, deltaS = self.getSolvationThermochemistryPierotti()
        return deltaH

    def getSolvationEnthalpyMintz(self):
        """Get the solvation enthalpy, in J/mol, from Mintz's expression for linear alkanes.
        
        See equation 15 in
        Enthalpy of Solvation Correlations For Gaseous Solutes
        Dissolved in Linear Alkanes (C5 - C16) Based on the Abraham Model
        C Mintz, K Burton, WE Acree Jr, MH Abraham
        http://dx.doi.org/10.1002/qsar.200730040
        """
        deltaH = -6.708 + 2.999*self.AbrahamE - 9.279*self.AbrahamL # kJ/mol
        return deltaH * 1000 # to get into J/mol
        
    def getSolvationEntropyMintz(self):
        """Get the solvation entropy, in J/mol, from Abraham's G and Mintz's H"""
        deltaH, deltaS = self.getSolvationThermochemistryMintz()
        return deltaS

    # Define which of the methods is returned as a property attribute
    SolvationFreeEnergy298 = property(getSolvationFreeEnergy298)
    SolvationThermochemistry = property(getSolvationThermochemistryMintz)
    SolvationEntropy = property(getSolvationEntropyMintz)
    SolvationEnthalpy = property(getSolvationEnthalpyMintz)
    
    
class PropertiesStore():
    """
    A class to store and evaluate Properties of all the Species in the system.
    """
    def __init__(self, resultsDir='RMG_results', speciesnames=None):
        self._specs_props = dict()
        self._speciesnames = speciesnames
        self.loadPropertiesFromFile(resultsDir)

    def __getattr__(self, property_name):
        """
        Get an array of the property. Length is Nspecies; 
        The order is the same as speciesnames array passed in at initialization, or 
        the order they are read from the file if no array is passed in.
        
        Uses self.getPropertyArray()
        
        >> s = SpeciesProperties()
        >> s.Radius
        
        Also works for things that are functions!
        >> s.PartitionCoefficientT(298)
        """
        try:
            return self.getPropertyArray(property_name, self._speciesnames)
        except (TypeError):
            try: 
                def fetcher_function(*args):
                    return self.getPropertyArray(property_name, self._speciesnames, *args)
                return fetcher_function
            except (KeyError, TypeError, ValueError), e:
                raise AttributeError
        except (KeyError, TypeError, ValueError), e:
            raise AttributeError
            
    def _getAttributeNames(self):
        """Get a list of extra (fake) attributes. Useful for ipython shell.""" 
        full_list = dir(PropertiesOfSpecies)
        return [a for a in full_list if not (a.startswith('__') or a.startswith('get')) ]
        
        
    def __getitem__(self,species_name):
        """Get PropertiesOfSpecies object of an individual species.
        
        Substitutes O2(1) for species without known properties (Ar and N2)
        then returns the corresponding PropertiesOfSpecies instance."""
        
        if species_name in ['Ar','N2']:
            # import pdb; pdb.set_trace()
            print "WARNING: Using 'O2(1)' properties for %s because I don't have %s values"%(species_name,species_name)
            species_name='O2(1)'
            
        spec_prop = self._specs_props[species_name]
        return spec_prop
        
    def loadPropertiesFromFile(self, resultsDir):
        """
        Load the species properties from the RMG_Solvation_Properties.txt file
        
        This will fill the self.properties dictionary, 
        with dictionaries of properties for each species.
        """
        if self._speciesnames is None:
            save_order = True # use the RMG_Solvation_Properties.txt file to determine species order 
            self._speciesnames = list()
        else:
            save_order = False # will leave the existing list speciesnames alone and use that for the order.
        import csv
        print "Loading species properties from",resultsDir
        propsfilename = os.path.join(resultsDir,'RMG_Solvation_Properties.txt')
        propsfile = open(propsfilename)
        reader = csv.DictReader(propsfile, dialect=csv.excel_tab)
        for spec_prop in reader:
            self._specs_props[spec_prop['ChemkinName']] = PropertiesOfSpecies(spec_prop)
            if save_order: self._speciesnames.append(spec_prop['ChemkinName'])
        propsfile.close()
        
    def getSpeciesProperty(self,species_name,property_name, *arguments):
        """
        Get the value of a property of a species. General method.
        
        Substitutes O2(1) for species without known properties (Ar, N2)
        then gets the attribute/property of the appropriate species.
        """
        if species_name in ['Ar','N2']:
            # import pdb; pdb.set_trace()
            # print "WARNING: Using 'O2(1)' properties for %s because I don't have %s values"%(species_name,species_name)
            species_name='O2(1)'
        try:
            spec_prop = self._specs_props[species_name]
        except KeyError:
            print "Don't have any properties for the species named '%s'."%species_name
            print "These are the species I have:",self._specs_props.keys()
            raise 
        #try:

        value = getattr(spec_prop,property_name)
        if isinstance(value,types.MethodType): # it's a function (method) not a value. call it and get the value
            value = value(*arguments)
        #except AttributeError:
            # print "Don't have the property '%s' for species '%s'."%(property_name,species_name)
            #raise
        return value
    
    def getPropertyArray(self,property_name, speciesnames, *arguments):
        """
        Get an array of the property values ordered by speciesnames. 
        
        Length is Nspecies; order is determined by speciesnames list.
        Probably quite slow, so wise to store the result."""
        values = numpy.zeros(len(speciesnames))
        
        for species_index, species_name in enumerate(speciesnames):
            values[species_index] = self.getSpeciesProperty(species_name,property_name, *arguments)
        return values
        


if __name__ == "__main__":
    import sys
    
    # use different chemistry mechanism if specified on the command line
    if len(sys.argv)>1:
        props = PropertiesStore(resultsDir = sys.argv[1])
    else:
        props = PropertiesStore()

