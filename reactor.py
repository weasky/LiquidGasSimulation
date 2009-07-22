"""
Get species and reactions info from .cti mechanism using ctml_writer
only NASA thermo is supported now
"""

import ctml_writer as ctml
from numpy import array,zeros

class OutsideValidRangeError(Exception):
    pass

units = ctml.units
OneAtm = ctml.OneAtm

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
        """
        getGibbsFreeEnergy(T) returns FreeEnergyOverR at temperature T
        """
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
        """ getThermo(T) returns (HeatCapacityOverR, EnthalpyOverRT, 
        EntropyOverR) for a given temperature T. 
        Raises outsideValidRangeError exception if T is not within range
        of polynomial.
        """
        if not self.ValidTemperature(T): raise outsideValidRangeError
        if self._pref <> 1.0e5: 
            raise Exception("not sure what to do with customised \
                    standard state pressure")
        # http://www.me.berkeley.edu/gri-mech/data/nasa_plnm.html
        # Cp/R = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4
        # H/RT = a1 + a2 T /2 + a3 T^2 /3 + a4 T^3 /4 + a5 T^4 /5 + a6/T
        # S/R  = a1 lnT + a2 T + a3 T^2 /2 + a4 T^3 /3 + a5 T^4 /4 + a7
        c=self._coeffs
        HeatCapacityOverR=0.0
        EnthalpyOverRT=0.0
        EntropyOverR=0.0
        for i in range(5):
            HeatCapacityOverR+= c[i] * T**i
            EnthalpyOverRT+= c[i] * T**i / (i+1)
        EnthalpyOverRT+=c[5]/T
        EntropyOverR = c[0]*math.log(T) + c[1]*T + c[2]*T*T/2 + \
                       c[3]*T*T*T/3 + c[4]*T*T*T*T/4 + c[6]
        return(HeatCapacityOverR,EnthalpyOverRT,EntropyOverR)
        
    def getGibbsFreeEnergy(self,T):
        """
        getGibbsFreeEnergy(T) returns  FreeEnergyOverR at temperature T
        """
        (HeatCapacityOverR,EnthalpyOverRT,EntropyOverR)= \
                                            self.getThermo(T)
        FreeEnergyOverRT=EnthalpyOverRT-EntropyOverR
        return FreeEnergyOverRT*T



class reaction(ctml.reaction):
    """A chemical reaction."""

    def existSpecies(self):
        for s in self._r.keys():
            if s not in ctml._speciesnames: return False
        return True

    def getForwardRateCoefficient(self,T):
        """returns the forward rate coefficient at a given temperature 
        assuming kf = [A, n, E] and kf=A (T/1K)^n exp(-E/RT)"""
        if not self.existSpecies: 
            raise Exception("Unknow species")
        kf=self._kf
        forwardRateCoefficient=kf[0] 
        if kf[1]:
            forwardRateCoefficient *= T**kf[1]
        assert ctml._ue=='kcal/mol', \
                                'unit convrsion not yet implemented'
        R=0.0019872065 # kcal/mol/K
        forwardRateCoefficient*=math.exp(-kf[2]/(R*T))
        print "%s forward k = %s"%(self._e,forwardRateCoefficient)
        return forwardRateCoefficient
        
    def getReverseRateCoefficient(self,T):
        """returns the reverse rate coefficient at given temperature"""
        if not self.rev: return 0

        reverseRateCoefficient=self.getForwardRateCoefficient(T)/ \
                self.getEquilibriumConstant(T)
        print "%s reverse k = %s"%(self._e,reverseRateCoefficient)
        return reverseRateCoefficient
        
    def getForwardRate(self,T,concentrations):
        """returns the forward rate of progress, given a temperature and a dictionary of species concentrations"""
        forwardRate=self.getForwardRateCoefficient(T)
        for speciesName,order in self._rxnorder.items():
        #   print "forwardRate*=concentrations[%s]**order = %s**%g"%(speciesName,concentrations[speciesName],order)
            forwardRate=forwardRate*concentrations[speciesName]**order
        print "%s forward rate = %s"%(self._e,forwardRate.simplified)
        return forwardRate
    
    def getReverseRate(self,T,concentrations):
        """returns the reverse rate of progress, given a temperature and a dictionary of species concentrations"""
        if not self.rev: return 0. # reaction not reversible
        reverseRate=self.getReverseRateCoefficient(T)
        for speciesName,order in self._p.items():
            #print "[%s]**%d=%s**%d"%(speciesName,order,concentrations[speciesName],order)
            reverseRate=reverseRate*concentrations[speciesName]**order
        print "%s reverse rate = %s"%(self._e,reverseRate.simplified)
        return reverseRate
    
    def getNetRate(self,T,concentrations):
        """returns the net rate of progress, given a temperature and a dictionary of species concentrations"""
        return self.getForwardRate(T,concentrations).simplified - self.getReverseRate(T,concentrations).simplified

    def getDeltaG(self,T):
        """
        returns the change in gibbs free energy *over R* at a given T
        """
        deltaGOverR=0
        for speciesName,order in self._p.items(): # add products
            deltaGOverR += order* \
            getSpeciesByName(speciesName).getGibbsFreeEnergy(T)
        for speciesName,order in self._r.items(): # subtract reactants
            deltaGOverR -= order* \
            getSpeciesByName(speciesName).getGibbsFreeEnergy(T)     
        return deltaGOverR
        
    def getEquilibriumConstant(self,T):
        """
        returns the equilibrium constant at a given temperature:  
        Keq = exp(-DeltaG/RT) 
        """
        
        #if not sum(self._p.values())==sum(self._rxnorder.values()): print "equilibrium constant calculation currently assumes forward and reverse reactions have same reaction order"
        
        dGoverR=self.getDeltaG(T) # returns delta G over R
        Keq= math.exp(-dGoverR/T) 
        
        # Kc=Keq * _uConc ** self.getDeltaNu()
        # if self.getDeltaNu(): print "this is wrong. Kc should not equal %s"%Kc
        # the units are right, but it's not _uConc that you should be multiplying by
        
        return Kc
        
    
    def getReactantNu(self):
        """Returns the stoichiometry in the forwards direction.
        
        i.e. the number of reactant molecules.
        uses self._reactantNu to cache the answer"""
        if hasattr(self,'_reactantNu'): return self._reactantNu
        reactantNu=0
        for speciesName,order in self._r.items(): # add reactants
            reactantNu += order
        self._reactantNu=reactantNu
        return reactantNu
    def getProductNu(self):
        """Returns the stoichiometry in the reverse direction. 
        
        i.e. the number of product molecules.
        uses self._productNu to cache the answer"""
        if hasattr(self,'_productNu'): return self._productNu
        productNu=0
        for speciesName,order in self._p.items(): # add products
            productNu += order
        self._productNu=productNu
        return productNu
    def getDeltaNu(self):
        """Returns the change in stoichiometry of the reaction, delta Nu.
        
        deltaNu= productNu - reactantNu
        uses self._deltaNu to cache the answer"""
        if hasattr(self,'_deltaNu'): return self._deltaNu
        self._deltaNu=self.getProductNu()-self.getReactantNu()
        return self._deltaNu


def getSpeciesByName(name):
    """Select a species by its name."""
    for s in ctml._species:
        if s._name==name:
            return s
    return None

def getNetRatesOfCreation(T,concs):
    """Get the net rate of creation of all the species at a given T 
    and concentration.
    
    Returns a dictionary. Expects a dictionary.
    not sure if should move this into reactor class"""
    nrocs=dict.fromkeys(_speciesnames)
    for speciesName in _speciesnames:
        nrocs[speciesName]=pq.Quantity(0.0,'mol/m**3/s')
    for r in _reactions:
        rate=r.getNetRate(T,concs)
        for speciesName,order in r._p.items(): # create products
            nrocs[speciesName]+=rate
            print "rate of creation of %s += %s. Now nroc=%s"%(speciesName,rate.simplified*order,nrocs[speciesName])
        for speciesName,order in r._r.items(): # consume reactants
            nrocs[speciesName]-=rate
            print "rate of consumption of %s += %s. Now nroc=%s"%(speciesName,rate.simplified*order,nrocs[speciesName])
    return nrocs # nrocs is a dictionary

    
class Reactor:
    def __init__(self, filename='',volume=1.0,
            pressure=OneAtm,temperature=400):
        if not ctml._species:
            execfile(filename)
        # default values
        self.t_start = 0.
        self.volume = volume
        self.P = pressure
        self.T = temperature
        self.nSpecies = len(ctml._species)
        self.nReactions = len(ctml._reactions)
        self.molWeight = zeros(self.nSpecies)
        self.massDens = zeros(self.nSpecies)
        self.speciesIndex = {}

    def setVolFrac(self,volFraction={}):
        """ set initial volume fraction:
        setIniVolFrac({'udecd':0.19,'xxx':0.2})
        todo: now all species involved should be included in the dict,
        even the conc is zero.
        """
        specienames = volFraction.keys()
        if len(specienames) <> self.nSpecies:
            raise Exception("lack species in input")
        for n,k in enumerate(specienames):
            self.speciesIndex[k]=n
        self.volFrac = volFraction.values()
        self.volFrac = array(self.volFrac)

    def setMolWeight(self,molWeight={}):
        """ set Molecular weight of species:
        reactor.setMolWeight({'udecd':135,'xxx':214})
        """
        for k,v in molWeight.items():
            i = self.speciesIndex[k]
            self.molWeight[i] = v

    def setMassDensity(self,massDensity={}):
        """ set density of species:
        reactor.setDensity({'udecd':710,'xxx':214})
        """
        for k,v in massDensity.items():
            i = self.speciesIndex[k]
            self.massDens[i] = v

    def getConc(self):
        """ update the concentrations in the reactor """
        self.molDens = self.massDens / self.molWeight 
        self.concs = self.volFrac * self.molDens
        return self.concs

    def setInitialTime(self, T0):
        self.t_start = T0

    def setT(self,T):
        """ set temperature of the reactor """
        self.T = T

    def getRightSideOfODE(self):
        """Basically the same as getNetRatesOfCreation() but takes 
        an array and returns an array
        since the order should be the same, we don't check the
        species' order in the dict  
        """
        concsDict = dict(zip(self.speciesIndex, self.concs))
        nrocsDict = getNetRatesOfCreation(self.T,concsDict)
        nrocsArray = array(nrocsDict.values())
        return nrocsArray

    def advance(self, time):
        """Deprecated.
        Advance the state of the reactor in time from the current
        time to time 'time'. Note: this method is deprecated. See
        class ReactorNet."""
        from scipy.integrate import odeint
        c0 = self.concs
        self.concs=odeint(self.getRightSideOfODE,c0,time)

    def step(self, time):
        """Deprecated.
        Take one internal time step from the current time toward
        time 'time'. Note: this method is deprecated. See class
        ReactorNet.
        todo"""
        raise "use method step of class ReactorNet"        
        #return _cantera.reactor_step(self.__reactor_id, time)    
    
        
if __name__ == "__main__":
    
    import sys, os
    if len(sys.argv)>1:
        file = sys.argv[1]
    else:
        print "using default file, chem.cti"
        file = 'chem.cti' # default input file
    base = os.path.basename(file)
    root, ext = os.path.splitext(base)
    dataset(root)

    test = Reactor(filename = file)
    test.setVolFrac({'O2(1)':0.01,'n-undecane(2)':0.04,
        'n-tridecane(3)':0.19,'SPC(4)':0.11,'n-hexadecane(5)':0.25,
        'SPC(6)':0.12,'n-nonadecane(7)':0.18,'n-heneicosane(8)':0.10,
        })
