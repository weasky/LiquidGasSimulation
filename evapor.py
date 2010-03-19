import ctml_writer as ctml
from numpy import pi,exp, append, arange, array, zeros, linspace
from scipy.integrate import odeint
from scipy.optimize import fsolve
from tools import diffusivity_hco_in_air

from chemistry import ChemistrySolver
import copy
import time

units = ctml.units
OneAtm = ctml.OneAtm
R = 8.314472 # J/mol K

class FuelComponent():
    """copy from Richard's chem solver
    specific to the fuel surrogate model.
    don't get confused with the reaction products."""
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
        self.liquidMolarDensity=liquidMolarDensity # mol/m3
        self.initialConcentration=self.liquidMolarDensity*initialVolFraction
        
    def getPsat(self, Temperature):
        """Use Antoine Equation to get saturated vapor pressure at Temperature.
        
        Note the units in Antoine Eqn is mmHg and C
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

class LiquidFilmCell:
    """
    This is a class for the thin liquid film on the cylinder wall. The instance
    represents each cell after descritization. when initialize, nSpecies is
    the number of hydrocarbon species, N2 and O2 will be automatically added.
    """
    def __init__(self, fuel=[], solver=ChemistrySolver(), diameter=1, thickness=1, length=1,
                 T=400, P=OneAtm, reaction=True):
        #use nSpecies, the code will calculate molWeight based on nC, nH and nO
        #2 means [O2 N2]
        # 7 surrogate model
        fuel=[
            FuelComponent('n-C11(2)',  0.05,dict(C=11,H=24,O=0),dict(A=6.9722, B=1569.57, C=187.7  ),4945.0),
            FuelComponent('n-C13(3)',  0.19,dict(C=13,H=28,O=0),dict(A=7.00756,B=1690.67, C=174.22 ),4182.0),
            FuelComponent('Mnphtln(4)',0.11,dict(C=11,H=10,O=0),dict(A=7.03592,B=1826.948,C=195.002),3.7e3 ),
            FuelComponent('n-C16(5)',  0.25,dict(C=16,H=34,O=0),dict(A=7.02867,B=1830.51, C=154.45 ),3415.0),
            FuelComponent('C10bnzn(6)',0.12,dict(C=16,H=26,O=0),dict(A=7.8148, B=2396.8,  C=199.5736),2.6e3),
            FuelComponent('n-C19(7)',  0.18,dict(C=19,H=40,O=0),dict(A=7.0153, B=1932.8,  C=137.6  ),2889.0),
            FuelComponent('n-C21(8)',  0.10,dict(C=21,H=44,O=0),dict(A=7.0842, B=2054,    C=120.1  ),2729.0)
        ]
        self.nSpecies = solver.Nspecies
        self.speciesnames = solver.speciesnames
        self.molWeight = solver.properties.MolecularWeight/1000. #to kg/m^3
        #evapFlux out
        self.evapFlux = zeros(self.nSpecies)
        # area and vol
        self.dia = diameter
        self.thickness = thickness
        self.initial_thickness = copy.deepcopy(thickness) # in case it's not just a float!
        self.len = length
        self.vol = pi * self.dia * self.len * self.thickness
        self.area = pi * (self.dia - 2 * self.thickness) * self.len

        self.T = T
        self.P = P
        self.reaction = reaction
        self.Psat = zeros(self.nSpecies)
        self.massDens = zeros(self.nSpecies)
        self.molDens = zeros(self.nSpecies)
        self.massFrac = zeros(self.nSpecies)
        self.molFrac = zeros(self.nSpecies)
        self.volFrac = zeros(self.nSpecies)
        self.concs = zeros(self.nSpecies)
        self.nC = solver.properties.nC
        self.nH = solver.properties.nH
        self.nO = solver.properties.nO
        self.Dvi = zeros(self.nSpecies)
        self.Dvi = diffusivity_hco_in_air(T=self.T,p=self.P*1.e-5,nC=self.nC,
                                          nH=self.nH,nO=self.nO)
        # process from fuel input
        concs_dict=dict.fromkeys(solver.speciesnames)
        volFrac_dict=dict.fromkeys(solver.speciesnames)
        antoine_dict=dict.fromkeys(solver.speciesnames)
        molDens_dict = dict(zip(solver.speciesnames,1./solver.properties.MolarVolume))
        for speciesName in solver.speciesnames:
            concs_dict[speciesName]= 0.0
            volFrac_dict[speciesName] = 0.0
            antoine_dict[speciesName] = 0.0
        for component in fuel:
            concs_dict[component.name]=component.initialConcentration # mol/m3
            volFrac_dict[component.name] = component.initialVolFraction
            antoine_dict[component.name] = component.Antoine
            molDens_dict[component.name] = component.liquidMolarDensity
        self.concs = array([concs_dict[s] for s in self.speciesnames])
        self.volFrac = array([volFrac_dict[s] for s in self.speciesnames])
        self.molDens = array([molDens_dict[s] for s in self.speciesnames])
        #updating
        self.massDens = self.molDens * self.molWeight
        molFrac =self.volFrac * self.molDens
        self.molFrac = molFrac / sum(molFrac)
        massFrac = self.volFrac * self.massDens
        self.massFrac = massFrac / sum(massFrac)
        self.vaporConcs = self.concs/solver.properties.PartitionCoefficient
        self.Psat = sum(self.concs) * R * self.T/solver.properties.PartitionCoefficient
        # correct Psat for fuel components (uses Antoine equation to correct for T)
        for component in fuel:
            species_index = self.speciesnames.index(component.name)
            self.Psat[species_index] = component.getPsat(self.T)
        self.vaporConcs = self.Psat * self.molFrac / R / self.T
        self.vaporMassDens = self.vaporConcs * self.molWeight
        self.vaporMassDens[:3] = 0.
        #air
        self.airMolFrac = zeros(2)
        self.airP = zeros(2)
        self.airMolWeight = array([32.0, 28.0134]) / 1000.
        self.airMassDens = array([1.429, 1.251])
        #tmp, may delete in future
        self.solver = solver
        self.solver.setTemperature(self.T)
        self.solver.setConcentrations(self.concs)
        #get air pressure and concs
        self.update()
        self.dryTime = self.getDryTime()

    def getDryTime(self):
        """
        Find the drying time for this film by solving the equation. however, this time may
        cause physical problems. It would be safer to get a 1% of initial film thickness.
        0.001 is an initial guess.
        """
        drytime = fsolve(self.getThickness,0.001)
        return drytime

    def getThickness(self,time):
        """
        Get the thickness after a given time, by advancing (a copy of) the simulation.
        """
        film = copy.deepcopy(self)
        film.advance(array([0,time]))
        return film.thickness 
        
    def setEvapFlux(self, evapFlux):
        """Set the flux of species evaporating."""
        assert evapFlux.size==self.nSpecies, "Expected an array of size self.nSpecies=%d"%self.nSpecies
        self.evapFlux = evapFlux
        
    def setCHO(self, nC, nH, nO):
        """ set the number of atoms for each species"""
        if(len(nC) != (self.nSpecies)):
            print "The length of nC is not correct\n"
            return
        self.nC = array(nC)
        self.nH = array(nH)
        self.nO = array(nO)
        molWeight = 12.011 * self.nC + 1.008 * self.nH + 16 * self.nO
        self.molWeight = molWeight / 1000 # kg/mol
        self.Dvi = diffusivity_hco_in_air(T=473, p=self.P / 10. ** 5,
                             nC=self.nC, nH=self.nH, nO=self.nO)

    def massDens(self, massDens):
        """ set density of the system"""
        self.massDens = array(massDens)
        self.molDens = self.massDens / self.molWeight

    def setMolDens(self, molDens):
        """ set mole density of the system"""
        self.molDens = array(molDens)
        self.massDens = self.molDens * self.molWeight

    def setVolFrac(self, volFrac):
        self.volFrac = volFrac
        self.concs = self.volFrac * self.molDens
        """ mole fraction"""
        molFrac = volFrac * self.molDens
        a = sum(molFrac)
        self.molFrac = molFrac / a
        """ mass fraction"""
        massFrac = volFrac * self.massDens
        a = sum(massFrac)
        self.massFrac = massFrac / a

    def update(self):
        self.vol = pi * self.dia * self.len * self.thickness
        self.area = pi * (self.dia - 2 * self.thickness) * self.len
        # air partial pressure
        tmp = self.P-sum(self.Psat * self.molFrac)
        # O2 vol frac 20.9% of air
        self.airP[0] = 0.209 * tmp
        self.airP[1] = tmp-self.airP[0]
        # air mole fraction, O2 and N2
        self.airMolFrac[0] = 19.71 * 10 ** -4 * self.airP[0] * 10 ** -5
        tmp_benzene=exp(-6.05445-4.95673/(self.T/100))
        tmp_decane=exp(-6.8288+0.3404/(self.T/100))
        tmp_nitrogen = 0.8*tmp_decane+0.2*tmp_benzene
        self.airMolFrac[1] = tmp_nitrogen * self.airP[1] * 10 ** -5
        # rescale the hydrocarbons' mole fraction
        tot = sum(self.airMolFrac)
        self.molFrac = self.molFrac * (1-tot)
        # rescale the vol frac, concentration, mass frac etc.
        # for now, keep those things unchanged
        # hack, should use dict in future
        self.concs[1] = sum(self.concs)*self.airMolFrac[1]
        self.concs[2] = sum(self.concs)*self.airMolFrac[0]


    def getVaporDens(self):
        Pi = self.Psat * self.molFrac
        Ri = R / self.molWeight
        rhovi = Pi / Ri / self.T #kg/m3
        return rhovi

    def vaporDiff(self, Lv=1):
        """
        The surface flux at the interface. A simple approximation is
        Q_i = Dv_i rhov_i/Lv
        return an array
        """
        self.vaporConcs = self.Psat * self.molFrac / R / self.T
        self.vaporMassDens = self.vaporConcs * self.molWeight
        self.vaporMassDens[:3] = 0.
        rhovi = self.vaporMassDens
        if (sum(self.evapFlux) == 0):
            Qi = self.Dvi * rhovi / Lv
            Qi = Qi * self.dia / 4. / self.len
        elif (sum(self.evapFlux) > 0):
            Qi = self.Dvi * self.evapFlux
        else:
            #print 'evaporation flux <0 ,something wrong'
            Qi = self.Dvi * self.evapFlux
        
        return Qi

    def rightSideofODE(self, Y, t):
        """
        drho/dt=A_l/V_l Qi - A_l/V_l (rhoi - sum(Qi)/sum(rhoi)) + source term
        y is mass fractional density
        """
        # evaporation term
        massFracDens = Y[:-1]
        molFracDens = massFracDens / self.molWeight
        self.massFrac = massFracDens / sum(massFracDens)
        self.molFrac = molFracDens / sum(molFracDens)
        ratio = self.area / self.vol
        Qi = self.vaporDiff(Lv=self.dia)
        Q = sum(Qi)
        #reaction source term, turn the mole to mass frac dens
        reactConcs = self.solver.getRightSideOfODE(molFracDens,t)*self.molWeight if self.reaction \
            else zeros(self.nSpecies)
   
        dhdt = -1. / sum(massFracDens) * Q
        drhodt = -ratio * Qi - ratio * massFracDens * dhdt + reactConcs
        drhodt = append(drhodt, dhdt)
        self.update()
        return drhodt

    def advance(self, t, plotresult=False):
        y0 = self.concs * self.molWeight
        y0 = append(y0, self.thickness)
        yt = odeint(self.rightSideofODE, y0, t)
        if(plotresult):
            import matplotlib.pyplot as plt
            plt.figure()
            plt.semilogy(t, yt)
            plt.legend(self.speciesnames)
            plt.show()
        self.thickness = yt[-1][-1]
        ytt = yt[-1][:-1]
        #        for iii in range(len(ytt)):
        #            if ytt[iii]<0:
        #                ytt[iii]=0.
        molFracDens = ytt / self.molWeight
        self.concs = molFracDens
        self.molFrac = molFracDens / sum(molFracDens)
        self.massFrac = ytt / sum(ytt)

if __name__ == "__main__":
    dia = 0.14E-3
    L = 0.5E-3
    initial_film_thickness = 3E-6

    diesel = LiquidFilmCell(T=473, diameter=dia, length=L, thickness=initial_film_thickness, reaction=False)
    # diesel.setCHO(nC=[11, 13, 11, 25, 16, 18, 10],
        #           nH=[24, 28, 10, 34, 26, 40, 44],
        #           nO=[0, 0, 0, 0, 0, 0, 0])
    print 'diesel components mol weight is', diesel.molWeight #g/mol
    # diesel.setMolDens(molDens=[4945.0, 4182.0, 3700.0,
        #               3415.0, 2600.0, 2889., 2729.]) #mol/m3
    print 'diesel components mol density is', diesel.molDens #
    print 'diesel components mass density is', diesel.massDens #kg/m3
    print 'the mol fraction is ', diesel.molFrac
    print 'the mass fraction is ', diesel.massFrac
    print 'the concentrations are ', diesel.concs
    print 'the total vapor pressure using K is ',sum(diesel.concs/diesel.solver.properties.PartitionCoefficient)*R*diesel.T
    print 'the saturated vapor pressure is', diesel.Psat
    print 'the total vapor pressure using Antoine is ',sum(diesel.Psat*diesel.molFrac)
    print 'the O2 and N2 partial pressure is ', diesel.airP[0], diesel.airP[1]
    print 'the O2 and n2 mole fraction are ', diesel.airMolFrac[0], diesel.airMolFrac[1]
        
    print 'the vapor densities are ', diesel.vaporMassDens
    qi = diesel.vaporDiff(Lv=dia)
    print 'the mass flux out of the interface ',qi
    print 'the initial h is', diesel.thickness
        # chem solver exp

    print 'initial concentrations are',diesel.concs
    print 'start evaporating without reaction'
    timesteps=linspace(0,0.5,501)
    diesel.advance(timesteps,plotresult=True)
    print 'the concentrations are ', diesel.concs
    print 'the vapor densities are ', diesel.getVaporDens()
    print 'the new h is', diesel.thickness
    print '%f percent film left', diesel.thickness / initial_film_thickness

    print 'start evaporating with reaction'
    diesel2 = LiquidFilmCell(T=473, diameter=dia, length=L, thickness=initial_film_thickness)
    timesteps=linspace(0,0.5,501)
    diesel2.advance(timesteps,plotresult=True)
    print 'the concentrations are ', diesel2.concs
    print 'the vapor densities are ', diesel2.getVaporDens()
    print 'the new h is', diesel2.thickness
    print '%f percent film left', diesel2.thickness / initial_film_thickness
