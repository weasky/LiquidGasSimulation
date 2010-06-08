import ctml_writer as ctml
from numpy import pi,exp, append, arange, array, zeros, linspace
import numpy
from scipy.integrate import odeint
from scipy.optimize import fsolve

from chemistry import ChemistrySolver
import copy
import time
import pylab

from PyDAS import dassl

from properties import PropertiesOfSpecies, PropertiesStore

units = ctml.units
OneAtm = ctml.OneAtm
R = 8.314472 # J/mol K

class FuelComponent():
    """
    A component of the surrogate diesel fuel.
    
    Don't get confused with the reaction products.
    """
    def __str__(self):
        return "<Species %s>"%(self.name)
    def __init__(self, name="",
            initialVolFraction=0,
            composition=dict(C=0,H=0,O=0), 
            Antoine=dict(A=0,B=0,C=0),
            liquidMolarDensity=3500 ):
        self.name = name
        self.initialVolFraction = initialVolFraction
        self.composition = composition
        self.Antoine = Antoine
        self.liquidMolarDensity = liquidMolarDensity # mol/m3
        self.initialConcentration = self.liquidMolarDensity*initialVolFraction
        
    def getPsat(self, Temperature):
        """
        Use Antoine Equation to get saturated vapor pressure at Temperature.

        Temperature is specified in K.
        Pressure is returned in Pa.
              
        Note, however, that the units in the Antoine parameters are mmHg and C
        P = 10^(A-B/(C+T))
        """
        A = self.Antoine['A']
        B = self.Antoine['B']
        C = self.Antoine['C']
        # transfer from mmHg to Pa
        # A = A0 + log10(101325/760)
        A = A + 2.124903
        # transfer from C to K
        C = C-273.15
        # Antoine's Equation
        return 10**(A-B / (C + Temperature))

class Phase(object):
    """A base class for a phase contaning some amounts of some species."""
    
    def __init__(self, properties, amounts=None):
        """Store the properties store instance, and cache some useful properties."""
        self.properties = properties
        # save these to avoid recomputing them
        self.molar_masses = self.properties.MolecularWeight/1000. #to kg 
        self.molar_volumes = self.properties.MolarVolume
        self.molar_densities = 1 / self.molar_volumes
        self.mass_densities = self.molar_densities * self.molar_masses 

        if amounts is None:
            self.amounts = zeros(properties.nSpecies, numpy.float64)
        else:
            self.amounts = amounts
            
    def get_total_amount(self):
        """The total amount of stuff in the phase, in mol."""
        return sum(self.amounts)
    total_amount = property(get_total_amount)
    
    def get_mole_fractions(self):
        """The amount (mole) fractions in the phase. Sum to 1."""
        return self.amounts / sum(self.amounts)
    mole_fractions = property(get_mole_fractions)

    def get_total_volume(self):
        """The total volume of the phase, in m3."""
        return sum(self.amounts * self.molar_volumes)
    total_volume = property(get_total_volume)
    
    def get_concentrations(self):
        """The concentrations of each species, in mol/m3."""
        if self.amounts.any():
            return self.amounts / self.get_total_volume()
        else: # they're all zero
            return self.amounts * 0
    concentrations = property(get_concentrations)
        
    def get_total_concentration(self):
        """The total molar concentration of the phase, in mol/m3."""
        return sum(self.amounts) / self.get_total_volume()
    total_concentration = property(get_total_concentration)
    
    
class DepositPhase(Phase):
    """
    A phase-separated deposit-forming phase inside the liquid film.
    """
    
    def __init__(self, properties, temperature, amounts=None):
        """
        * properties is a PropertiesStore instance.
        * temperature is in Kelvin.
        * amounts is an array, storing the number of moles of each species 
          inside this phase; zeros if unspecified.
        """
        Phase.__init__(self, properties, amounts)
        self.T = temperature
        self.DETERGENCY = None

        """
        Liquid-liquid mass transfer coefficients for stirred systems 
        seem to be of the order of magnitude 1e-5 m/s (eg. http://www.kirj.ee/public/va_ke/k50-1-3.pdf)
        which seems to be roughly what the diffusivities in air are (in m2/s),
        so we shall just use these for now (they're probably correlated).
        """
        self.mass_transfer_coefficients = self.properties.DiffusivityInAir(self.T, 1e5) # m2/s, though we're using as m/s
        self.mass_transfer_coefficients *= 1e-3 # slow it down for testing
    
    def get_total_volume(self):
        """The total volume of the phase, in m3.
        
        This method over-rides the parent class one which is based on molar volumes.
        """
        return self.initial_volume
    total_volume = property(get_total_volume)
    
    def get_diesel_equilibrium_concentrations(self):
        """Solution-phase concentrations at interface at equilibrium, in mol/m3.
        
        Currently this is not temperature-dependent.
        """
        DETERGENCY = self.DETERGENCY if self.DETERGENCY else 1 
        
        K = self.properties.DepositPartitionCoefficient298
        K = numpy.exp( numpy.log(K) / DETERGENCY )
        diesel_concentrations = self.concentrations * K
        return diesel_concentrations
    diesel_equilibrium_concentrations = property(get_diesel_equilibrium_concentrations)

    def add_detergent(self,capacity=2):
        """
        Simulate the one-parameter(capacity) detergent effect.

        From Richard's code, it seems
        K^C = K_0 (K_0 is partitioncoef without detergent, C is capacity)

        that means detergent no effect if C=1 and Large capacity will give
        smaller K if C>1
        """
        self.DETERGENCY = capacity
        

    def fluxes_in(self, outside_concentrations):
        """The fluxes of species into the deposit phase from diesel at outside_concentration, in mol/m2/s."""
        driving_forces = outside_concentrations - self.diesel_equilibrium_concentrations
        fluxes = driving_forces * self.mass_transfer_coefficients
        return fluxes


    def equilibrate(self, outside_concentrations):
        """
        Bring into equilibrium with the outside.
        
        Updates self.amounts and returns the amounts on the outside after equilibration.
        
        THIS FUNCTION IS WRONG
        """
        pass
     #   total_amounts = outside_amounts + self.amounts
     #   # split is the ratio between outside and inside the deposit phase
     #   split = self.properties.DepositPartitionCoefficient298
     #   
     #   # outside / self = split
     #   self.amounts = total_amounts / (split + 1)
     #   outside_amounts = self.amounts * split
     #   return outside_amounts



class LiquidFilmCell(dassl.DASSL, Phase):
    """
    This is a class for the thin liquid film on the cylinder wall. The instance
    represents each cell after discretization. When initializing, nSpecies is
    the number of hydrocarbon species, N2 and O2 will be automatically added.
    
    Temperature is constant throughout.
    """
    
    def __init__(self, chem_solver=None,
                 diameter=1, thickness=1, length=1,
                 T=400, P=OneAtm,
                 EVAPORATION=True,
                 CHEMICAL_REACTION=True,
                 PHASE_SEPARATION=True,
                 resultsDir='RMG_results'):
        # Initialize the DASSL solver
        dassl.DASSL.__init__(self)  

        
        # store parameters
        self.T = T
        self.P = P
        
        self.EVAPORATION = EVAPORATION
        self.CHEMICAL_REACTION = CHEMICAL_REACTION
        self.PHASE_SEPARATION = PHASE_SEPARATION
        
        # the vector 'y' solved by dassl is the amounts, divided by SCALE
        self.SCALE = 1e-12  # so dassl is tracking some numbers larger than the actual amounts 
        
        # Set up and store chemistry solver
        if chem_solver is None: chem_solver = ChemistrySolver(resultsDir=resultsDir)
        self.chem_solver = chem_solver
        self.nSpecies = chem_solver.Nspecies
        self.speciesnames = chem_solver.speciesnames
        assert len(self.speciesnames) == self.nSpecies        

        # Create properties store.
        # Must pass it list of speciesnames so that the arrays returned are in the right size and order
        self.properties = PropertiesStore(resultsDir=resultsDir, speciesnames=self.speciesnames)
        
        # 7 surrogate model
        if all(name in self.speciesnames for name in 'O2(1)  n-C11(2)  n-C13(3) \
                  Mnphtln(4)  n-C16(5)  C10bnzn(6)  n-C19(7)  n-C21(8)'.split() ):
            # We're using Richard's new species names
            fuel=[
                FuelComponent('n-C11(2)',  0.05,dict(C=11,H=24,O=0),dict(A=6.9722, B=1569.57, C=187.7  ),4945.0),
                FuelComponent('n-C13(3)',  0.19,dict(C=13,H=28,O=0),dict(A=7.00756,B=1690.67, C=174.22 ),4182.0),
                FuelComponent('Mnphtln(4)',0.11,dict(C=11,H=10,O=0),dict(A=7.03592,B=1826.948,C=195.002),3.7e3 ),
                FuelComponent('n-C16(5)',  0.25,dict(C=16,H=34,O=0),dict(A=7.02867,B=1830.51, C=154.45 ),3415.0),
                FuelComponent('C10bnzn(6)',0.12,dict(C=16,H=26,O=0),dict(A=7.8148, B=2396.8,  C=199.5736),2.6e3),
                FuelComponent('n-C19(7)',  0.18,dict(C=19,H=40,O=0),dict(A=7.0153, B=1932.8,  C=137.6  ),2889.0),
                FuelComponent('n-C21(8)',  0.10,dict(C=21,H=44,O=0),dict(A=7.0842, B=2054,    C=120.1  ),2729.0)
            ]
            self._oxygen_name = "O2(1)"
        elif all(name in self.speciesnames for name in "SPC(1)  SPC(2)  SPC(3) \
                              SPC(4)  SPC(5)  O2(6)  SPC(7)  SPC(8)".split() ):
            # We're using the old names (and had better hope they're in the right order!)
            fuel=[                                                                                                  
                FuelComponent('SPC(2)',  0.05,dict(C=11,H=24,O=0),dict(A=6.9722, B=1569.57, C=187.7  ),4945.0), # n-C11(2)
                FuelComponent('SPC(3)',  0.19,dict(C=13,H=28,O=0),dict(A=7.00756,B=1690.67, C=174.22 ),4182.0), # n-C13(3) 
                FuelComponent('SPC(4)',  0.11,dict(C=11,H=10,O=0),dict(A=7.03592,B=1826.948,C=195.002),3.7e3 ), # Mnphtln(4)
                FuelComponent('SPC(5)',  0.25,dict(C=16,H=34,O=0),dict(A=7.02867,B=1830.51, C=154.45 ),3415.0), # n-C16(5) 
                FuelComponent('SPC(1)',  0.12,dict(C=16,H=26,O=0),dict(A=7.8148, B=2396.8,  C=199.5736),2.6e3), # C10bnzn(6)
                FuelComponent('SPC(7)',  0.18,dict(C=19,H=40,O=0),dict(A=7.0153, B=1932.8,  C=137.6  ),2889.0), # n-C19(7) 
                FuelComponent('SPC(8)',  0.10,dict(C=21,H=44,O=0),dict(A=7.0842, B=2054,    C=120.1  ),2729.0)  # n-C21(8)
            ]
            self._oxygen_name = "O2(6)"
        else:
            raise Exception("Don't recognize species names used for fuel components.")

        # Initialize the Phase parent class; sets up 'amounts' variable and caches some properties
        Phase.__init__(self, self.properties)

        # Deposit
        self.deposit = DepositPhase(self.properties, self.T)
        
        # store some other unchanging properties here, to save time looking them up or calculating them repeatedly
        self.Dvi = self.properties.DiffusivityInAir(self.T, self.P) # m2/s


        # save some special species indices
        self._oxygen_index = self.speciesnames.index(self._oxygen_name)
        self._nitrogen_index = self.speciesnames.index('N2')
        
        # area and volume of liquid film cell
        self.diameter = diameter
        self.initial_thickness = copy.deepcopy(thickness) # in case it's not just a float!
        self.length = length
        self.initial_volume = pi * diameter * length * thickness
        
        # deposit volume (held constant and used for concentrations!)
        self.deposit.initial_volume = self.initial_volume * 0.01
        
        # Set the fuel component initial amounts
        # this enables us to find the total concentration (at the start) necessary for Psat calculations
        for component in fuel:
            species_index = self.speciesnames.index(component.name)
            self.amounts[species_index] = ( self.initial_volume * 
                                            component.initialVolFraction * 
                                            component.liquidMolarDensity )
        
        # Set the Psat values
        """
        Psat is not really a saturated vapor pressure, 
        but is the x=1 extrapolation of the partial pressure vs. x 
        line based on the partition coefficient from the Abraham model.
        
        In the Henry's Law model used later to get vapor phase concentration,
         Cvi = xsi * Psati / RT            (1)
        The vapor/solution partition coefficient Kvsi is defined as
         Kvsi = Csi / Cvi                  (2)
        Substituting (2) into (1), along with xi = Csi/Cstotal gives
         'Psati' = RT Cstotal / Kvsi
         
         We can only calculate this once we know the total concentration.
         We will assume it doesn't change significantly from the start, but
         we have to wait until the intial fuel concentrations have been set
         before we can evaluate these terms.
        """
        self.Psats = self.get_total_concentration() * R * self.T / self.properties.PartitionCoefficientT(self.T)
        
        # Update fuel component parameters with revised values.
        # (we have better information for the fuel components than for other species)
        for component in fuel:
            species_index = self.speciesnames.index(component.name)
            self.molar_volumes[species_index] = 1 / component.liquidMolarDensity
            # Psat for fuel components uses Antoine equation
            self.Psats[species_index] = component.getPsat(self.T)
        
        
        # Some species we don't want "evaporating" so set their Psat values to zero.
        for species_name in ['Ar', 'N2', self._oxygen_name]:
            species_index = self.speciesnames.index(species_name)
            self.Psats[species_index] = 0
        
       ## put some diesel in the deposit phase to get it started!
       #DEPOSIT_START_FRACTION = 0.01
       #self.deposit.amounts = self.amounts * DEPOSIT_START_FRACTION
       #self.amounts = self.amounts * (1 - DEPOSIT_START_FRACTION)
        
        # Set up the chemistry solver
        self.chem_solver.setTemperature(self.T)
        self.chem_solver.setConcentrations(self.concentrations)
        # get amount of air dissolved in the diesel
        self.update_oxygen_nitrogen()
        # find, through trial and error, how long it takes to dry out
        # self.dryTime = self.getDryTime()
        
    ############
    # Attribute-style parameter functions:
    ############
 
    def get_vapor_concentrations(self):
        """Vapor-phase concentrations at interface, in mol/m3."""
        vaporConcs = self.Psats * self.mole_fractions / R / self.T
        # an alternative would be to get them from 
        # vaporConcs = self.concs / self.properties.PartitionCoefficientT(self.T)
        return vaporConcs
    vapor_concentrations = property(get_vapor_concentrations)
    
    def get_vapor_partial_pressures(self):
        """Vapor-phase partial pressures at interface, in Pa."""
        partial_pressures = self.Psats * self.mole_fractions
        return partial_pressures
    vapor_partial_pressures = property(get_vapor_partial_pressures)
    
    def get_thickness(self):
        """Thickness of the film, in m."""
        thickness = self.total_volume / pi / self.diameter / self.length
        return thickness
    thickness = property(get_thickness)

    def get_area(self):
        """Area of the interface with the vapor phase, in m^2.
        
        Takes into account the thickness of the film, which is probably negligible and unnecessary.
        """
        return pi * (self.diameter - 2 * self.thickness) * self.length
    area = property(get_area)
    
    #############
    ## update functions
    #############

    def update_oxygen_nitrogen(self):
        """
        Update the amount of oxygen and nitrogen dissolved in the deposit.
        
        This uses some relationship between mole fraction and partial pressure,
        that Yinchun has found somewhere. Please could you elaborate?
        """
        air_partial_pressure = self.P - sum(self.vapor_partial_pressures)
        assert air_partial_pressure > 0
        oxygen_partial_pressure = 0.209 * air_partial_pressure
        nitrogen_partial_pressure = air_partial_pressure - oxygen_partial_pressure
        # air mole fraction, O2 and N2. 
        # Not sure where these correlations came from. Perhaps Yinchun can help..?
        oxygen_mole_fraction = 19.71E-4 * oxygen_partial_pressure * 1E-5
        
        tmp_benzene=exp(-6.05445-4.95673/(self.T/100))
        tmp_decane=exp(-6.8288+0.3404/(self.T/100))
        tmp_nitrogen = 0.8*tmp_decane+0.2*tmp_benzene
        nitrogen_mole_fraction = tmp_nitrogen * nitrogen_partial_pressure * 1E-5
        
        self.amounts[self._oxygen_index] = self.total_amount * oxygen_mole_fraction
        self.amounts[self._nitrogen_index] = self.total_amount * nitrogen_mole_fraction
    
    def get_oxygen_dissolution_rate(self):
        """The rate at which oxygen is dissolving into the diesel, in mol/s.
        
        We asssume the equilibrium mole fraction is given by Yinchun's correlation,
        and the rate at which this equilbrium is reached has a first order 
        decay constant arbitrarily chosen, eg. k=1e3 s^-1.
        """
        
        # from Yinchun's email dated 8 Feb 2010
        # Mole fractions in the hydrocarbon at 1 atm partial pressure,
        # multiplied by their partial pressures in pure air 
        # (not allowing for the hydrocarbons in the air)
        Xo2 = 21e-4 * 0.209 
       # Xn2 = 225e-4 * 0.79
        driving_force =  self.total_amount * Xo2 - self.amounts[self._oxygen_index]
        rate_constant = 1e2 # seconds^-1
        return driving_force * rate_constant
        
    ############
    ## legacy functions not needed
    ###########

    def getDryTime(self):
        """
        Find the drying time for this film 
        
        Works by repeatedly calling getThickness with different end times.
        """
        drytime = fsolve(self.getThickness,0.001)
        return drytime

    def getThickness(self,time):
        """
        Get the thickness after a given time, by advancing (a copy of) the simulation.
        """
        film = copy.deepcopy(self)
        film.advance(array([0,time], numpy.float64))
        return film.thickness 
        
    def get_vapor_mass_densities(self):
        """Vapor phase mass densites, in kg/m3"""
        Pi = self.Psats * self.mole_fractions
        Ri = R / self.molar_masses
        rhovi = Pi / Ri / self.T #kg/m3
        return rhovi
    
    #######
    ## simulation!!!
    #######
    
    def get_evaporative_flux(self, Lv=1):
        """
        The surface flux at the interface, in mol/m2/s.
        
        A simple approximation is
        Q_i = Dv_i rhov_i/Lv
        
        Where rhov_i is the molar density (concentration) in the vapor phase
        at the interface and Lv is a characteristic length.
        
        Returns an array.
        """        
        rhovi = self.get_vapor_concentrations()
        # Dvi : m2/s
        # rhovi : mol/m3
        # Lv: m
        # Qi: mol/m2/s
        Qi = self.Dvi * rhovi / Lv
        
        """
        Since the above flux is from the nozzle to the chamber, the A = \pi (D/2)^2
        but what we need if the flux from the fuel film, where A_film = \pi D L
        so we need a correction coeff D/(4L).

        Here we have an assumption that the rhovi doesn't change that much along
        the axis though
        """
        Qi = Qi * self.dia / 4. / self.len
        return Qi
    
    
    def residual(self, t, y, dydt):
        """The residual function, solved by dassl.
        
        Evaluate the residual function for this model, given the current value
        of the independent variable `t`, dependent variables `y`, and first
        derivatives `dydt`. Return a numpy array with the values of the residual
        function and an integer with status information (0 if okay, -2 to
        terminate).
        
        if:       dy/dt = g(t,y,dy/dt)
        then:  residual = g(t,t,dy/dt) - dy/dt
        
        In this case, y is the amounts of each species in the diesel phase, then
        the amounts of each species in the deposit phase.
        
        NOTE! This function returns a tuple: (residual, 0)
        """
        
        self.amounts = self.SCALE * y[:self.nSpecies] 
        self.deposit.amounts = self.SCALE * y[self.nSpecies:]
        
        # currently P_air < 0 because too many volatile species
        # and I think updating the amount in here is a bad idea anyway -
        # better to let the DAE solver know about the change.
        #self.update_oxygen_nitrogen() # changes the amounts of these
        
        
        if self.PHASE_SEPARATION:
            dNdt_into_deposit = self.area * self.deposit.fluxes_in(self.concentrations)
            dNdt_into_deposit[self._oxygen_index] = 0
            dNdt_into_deposit[self._nitrogen_index] = 0
            #print 'dNdt_into_deposit',repr(dNdt_into_deposit)
        else:
            dNdt_into_deposit = numpy.zeros_like(self.amounts)
        
        if self.CHEMICAL_REACTION:
            # chem_solver deals with concentrations; to get amounts, scale by total volume
            dNdt_reaction =  self.total_volume * self.chem_solver.getRightSideOfODE(self.concentrations)
            #print 'dNdt_reaction',dNdt_reaction
        else:
            dNdt_reaction = numpy.zeros_like(self.amounts)
            
        if self.EVAPORATION:
            # evaporative_flux is per surface area
            dNdt_evaporation = self.area * self.get_evaporative_flux(Lv=self.diameter)
            #print 'dNdt_evaporation',dNdt_evaporation
        else:
            dNdt_evaporation = numpy.zeros_like(self.amounts)
        
        
        dNdt_diesel = dNdt_reaction - dNdt_into_deposit - dNdt_evaporation
        
        # add on oxygen dissolution rate
        dNdt_diesel[self._oxygen_index] += self.get_oxygen_dissolution_rate()
        
        g = numpy.concatenate((dNdt_diesel,dNdt_into_deposit)) / self.SCALE
        #print 'g=',repr(g)
        #print "y=",repr(y)
        #print "dydt=",repr(dydt)
        residual = g - dydt
        #print "residual=",repr(residual)
        return (residual, 0)
        
    def initialize_solver(self, time=0, atol=1e-15, rtol=1e-8):
        """Initialize the DASSL solver."""
       # self.nonnegative = True
        y = numpy.concatenate((self.amounts,self.deposit.amounts)) / self.SCALE
        # if simple ODE not DAE then residual with dydt=0 is in fact dydt
        # the residual method returns (residual,0) so take [0] element to get residual
        dydt = self.residual(time, y, numpy.zeros_like(y))[0] 
        print "I have dydt0 = ",repr(dydt)
        dassl.DASSL.initialize(self, time, y0=y, dydt0=dydt, atol=atol, rtol=rtol)
        

    def rightSideofODE(self, Y, t):
        """
        drho/dt=A_l/V_l Qi - A_l/V_l (rhoi - sum(Qi)/sum(rhoi)) + source term
        y is mass density
        """
        # evaporation term
        massDens = Y[:-1]
        molDens = massDens / self.molar_masses
        massFrac = massDens / sum(massDens)
        molFrac = molDens / sum(molDens)
        self.amounts = molDens * self.volume
        ratio = self.area / self.volume
        Qi = self.get_evaporative_flux(Lv=self.diameter)  # I think this is broken 
        # since I changed evaporative_flux to actually be a flux (per unit area)
        Q = sum(Qi)
        
        #reaction source term, turn the mole to mass frac dens
        if self.reaction:
            reactConcs = self.chem_solver.getRightSideOfODE(molDens,t)*self.molar_masses
        else:
            reactConcs = zeros(self.nSpecies, numpy.float64)
   
        dhdt = -1. / sum(massDens) * Q
        drhodt = -ratio * Qi - ratio * massDens * dhdt + reactConcs
        drhodt = append(drhodt, dhdt)
        self.update_volume_area_thickness()
        self.update_oxygen_nitrogen()
        return drhodt

    #def advance(self, t, plotresult=False):
    #    y0 = self.concentrations * self.molar_masses
    #    y0 = append(y0, self.thickness)
    #    yt = odeint(self.rightSideofODE, y0, t)
    #    if(plotresult):
    #        import matplotlib.pyplot as plt
    #        plt.figure()
    #        plt.axes([0.1,0.1,0.6,0.85])
    #        plt.semilogy(t, yt)
    #        plt.ylabel('mass concentrations (kg/m3)')
    #        plt.xlabel('time(s)')
    #        #plt.legend(self.speciesnames)
    #        for i in range(len(self.speciesnames)):
    #            plt.annotate(self.speciesnames[i], (t[-1],yt[-1,i]), 
    #                xytext=(20,-5), textcoords='offset points', 
    #                arrowprops=dict(arrowstyle="-") )
    #        plt.show()
    #    #self.thickness = yt[-1][-1]
    #    ytt = yt[-1][:-1]
    #    #        for iii in range(len(ytt)):
    #    #            if ytt[iii]<0:
    #    #                ytt[iii]=0.
    #    molDens = ytt / self.molar_masses
    #    concentrations = molDens
    #    self.amounts = concentrations * self.volume
    #    

        

if __name__ == "__main__":
    dia = 0.14E-3
    L = 0.5E-3
    initial_film_thickness = 3E-6

    diesel = LiquidFilmCell(T=473, diameter=dia, length=L, thickness=initial_film_thickness,
                            EVAPORATION=False, CHEMICAL_REACTION=True, PHASE_SEPARATION=True,
                            resultsDir='RMG_results-amrit' )

    print 'diesel components molar mass is', diesel.molar_masses # kg/mol
    print 'diesel components molar density is', diesel.molar_densities # mol/m3
    print 'diesel components mass density is', diesel.mass_densities # kg/m3
    print 'the mol fraction is', diesel.mole_fractions
    print 'the concentrations are', diesel.concentrations
    print 'the total vapor pressure using K is',sum(diesel.concentrations/diesel.properties.PartitionCoefficientT(diesel.T))*R*diesel.T
    print 'the saturated vapor pressure is', diesel.Psats
    print 'the total vapor pressure using Antoine is ',sum(diesel.Psats*diesel.mole_fractions)
    print 'the total vapor pressure used in the model is ',sum(diesel.vapor_partial_pressures)
    
    print 'the vapor densities are ', diesel.get_vapor_mass_densities()
    qi = diesel.get_evaporative_flux(Lv=dia)
    print 'the mass flux out of the interface ',qi
    print 'the initial h is', diesel.thickness
    print 'initial concentrations are',diesel.concentrations
    
    # import pdb; pdb.set_trace(); # pause for debugging
    if True:
        print "Trying DASSL solver"
        diesel.initialize_solver()
        timesteps=linspace(0,1,1001)
        
        #check the residual works (although the initialise_solver above has just done so)
        #import pdb; pdb.set_trace()
        residual = diesel.residual(diesel.t, diesel.y, diesel.dydt)
        
        diesel_history_array = numpy.zeros((len(timesteps),diesel.nSpecies), numpy.float64)
        deposit_history_array = numpy.zeros((len(timesteps),diesel.nSpecies), numpy.float64)
        for step,time in enumerate(timesteps):
            if time>0 : diesel.advance(time)
            diesel_history_array[step] = diesel.y[:diesel.nSpecies] * diesel.SCALE
            deposit_history_array[step] = diesel.y[diesel.nSpecies:] * diesel.SCALE
            print diesel.t, diesel.y
    
    
    final_amounts = diesel_history_array[-1]
    pylab.figure(1)
    pylab.axes([0.1,0.1,0.71,0.85])
    pylab.plot(timesteps,diesel_history_array)
    for i in range(len(diesel.speciesnames)):
        pylab.annotate(diesel.speciesnames[i], (timesteps[-1]*0.999,final_amounts[i]), 
            xytext=(20,-5), textcoords='offset points', 
            arrowprops=dict(arrowstyle="-") )
    pylab.title('diesel amounts')
    pylab.show()
    
    final_amounts = deposit_history_array[-1]
    pylab.figure(2)
    pylab.axes([0.1,0.1,0.71,0.85])
    pylab.plot(timesteps,deposit_history_array)
    for i in range(len(diesel.speciesnames)):
        pylab.annotate(diesel.speciesnames[i], (timesteps[-1]*0.999,final_amounts[i]), 
            xytext=(20,-5), textcoords='offset points', 
            arrowprops=dict(arrowstyle="-") )
    pylab.title('deposit amounts')
    pylab.show()
        
        
    #else:
    #    print 'start evaporating without reaction'
    #    timesteps=linspace(0,0.5,501)
    #    diesel.advance(timesteps,plotresult=True)
    #    print 'the concentrations are ', diesel.concentrations
    #    print 'the vapor densities are ', diesel.get_vapor_mass_densities()
    #    print 'the new h is', diesel.thickness
    #    print 'fraction of film left', diesel.thickness / initial_film_thickness
    #    
    #    print 'start evaporating with reaction'
    #    diesel2 = LiquidFilmCell(T=473, diameter=dia, length=L, thickness=initial_film_thickness)
    #    timesteps=linspace(0,0.5,501)
    #    diesel2.advance(timesteps,plotresult=True)
    #    print 'the concentrations are ', diesel2.concentrations
    #    print 'the vapor densities are ', diesel2.get_vapor_mass_densities()
    #    print 'the new h is', diesel2.thickness
    #    print 'fraction of film left', diesel2.thickness / initial_film_thickness
        
