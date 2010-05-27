#
# load, store and calculate properties of species in the simulation
# 
# Richard West 2010

import sys, os
import math
import pylab, numpy
import re

    
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
    PartitionCoefficient = property(getPartitionCoefficient)
    
    
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

    def getSolvationThermochemistry(self):
        """
        This is the Java code from RMG:
        
        double r_solute=p_chemGraph.getRadius();
        double r_solvent; r_solvent=3.498*Math.pow(10,-10);// 3.311;  // Manually assigned solvent radius [=] meter Calculated using Connolly solvent excluded volume from Chem3dPro
        double r_cavity=r_solute+r_solvent;                    // Cavity radius [=] Angstrom
        double rho; rho=0.00309*Math.pow(10,30);  //0.00381;   // number density of solvent [=] molecules/Angstrom^3   Value here is for decane using density =0.73 g/cm3
        double parameter_y=4.1887902*rho*Math.pow(r_solvent, 3);  // Parameter y from Ashcraft Thesis Refer pg no. 60. (4/3)*pi*rho*r^3
        double parameter_ymod=parameter_y/(1-parameter_y);     // parameter_ymod= y/(1-y) Defined for convenience
        double R=8.314;                                        // Gas constant units J/mol K
        double T=298;                                          // Standard state temperature
        
        // Definitions of K0, K1 and K2 correspond to those for K0', K1' and K2' respectively from Ashcraft's Thesis (-d/dT of K0,K1,K2)
        double K0= -R*(-Math.log(1-parameter_y)+(4.5*parameter_ymod*parameter_ymod));
        double K1= (R*0.5/r_solvent)*((6*parameter_ymod)+(18*parameter_ymod*parameter_ymod));
        double K2= -(R*0.25/(r_solvent*r_solvent))*((12*parameter_ymod)+(18*parameter_ymod*parameter_ymod));
        
        // Basic definition of entropy change of solvation from Ashcfrat's Thesis
        double deltaS0;
        deltaS0=K0+(K1*r_cavity)+(K2*r_cavity*r_cavity);
        
        // Generation of Abraham Solute Parameters
        AbramData result_Abraham= new AbramData();
        result_Abraham= p_chemGraph.getAbramData();
        
        // Solute descriptors from the Abraham Model
        double S=result_Abraham.S;
        double B=result_Abraham.B;
        double E=result_Abraham.E;
        double L=result_Abraham.L;
        double A=result_Abraham.A;
        
        //Manually specified solvent descriptors 
        // (constants here are for dry decane, from from M.H. Abraham et al. / J. Chromatogr. A 1037 (2004) 29–47 )
        double c = 0.156;
        double s = 0;
        double b = 0; 
        double e = -0.143;
        double l = 0.989;
        double a = 0;
        
        // Dry octan-1-ol, from M.H. Abraham et al. / J. Chromatogr. A 1037 (2004) 29–47 
        // c=-0.120; e=-0.203; s=0.560; a=3.560; b=0.702; l=0.939
        
        double logK = c + s*S + b*B + e*E + l*L + a*A;    // Implementation of Abraham Model
        double deltaG0_decane = -8.314*298*logK;
        //       System.out.println("The free energy of solvation in decane at 298K w/o reference state corrections  = " + deltaG0_decane +" J/mol for " );
        
        // Calculation of enthalpy change of solvation using the data obtained above
        double deltaH0=deltaG0_decane+(T*deltaS0);
        deltaS0=deltaS0/4.18;   //unit conversion from J/mol to cal/mol
        deltaH0=deltaH0/4180;   //unit conversion from J/mol to kcal/mol
        """
        pass
            
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
        
        >> s = SpeciesProperties()
        >> s.Radius
        """
        try:
            return self.getPropertyArray(property_name, self._speciesnames)
        except KeyError:
            raise AttributeError
            
    def _getAttributeNames(self):
        """Get a list of extra (fake) attributes. Useful for ipython shell.""" 
        full_list = dir(PropertiesOfSpecies)
        return [a for a in full_list if not (a.startswith('__') or a.startswith('get')) ]
        
        
    def __getitem__(self,species_name):
        """Get properties of an individual species."""
        
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
        
    def getSpeciesProperty(self,species_name,property_name):
        """Get the value of a property of a species. General method."""
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
        #except AttributeError:
            # print "Don't have the property '%s' for species '%s'."%(property_name,species_name)
            #raise
        return value
    
    def getPropertyArray(self,property_name, speciesnames):
        """Get an array of the property. Length is Nspecies; order is same as chemistry model.
        
        Probably quite slow, so wise to store the result."""
        values = numpy.zeros(len(speciesnames))
        for species_index, species_name in enumerate(speciesnames):
            values[species_index] = self.getSpeciesProperty(species_name,property_name)
        return values


if __name__ == "__main__":
    import sys
    
    # use different chemistry mechanism if specified on the command line
    if len(sys.argv)>1:
        props = PropertiesStore(resultsDir = sys.argv[1])
    else:
        props = PropertiesStore()

