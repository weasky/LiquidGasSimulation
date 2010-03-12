from math import pi,exp

import ctml_writer as ctml
from numpy import append
from numpy import arange
from numpy import array
from numpy import zeros
from scipy.integrate import odeint
from tools import diffusivity_hco_in_air

units = ctml.units
OneAtm = ctml.OneAtm
R = 8.314472 # J/mol K


class LiquidFilmCell:
    """
    This is a class for the thin liquid film on the cylinder wall. The instance
    represents each cell after descritization. when initialize, nSpecies is
    the number of hydrocarbon species, N2 and O2 will be automatically added.
    """
    def __init__(self, nSpecies=1, diameter=1, thickness=1, length=1,
                 T=400, P=OneAtm):
        #use nSpecies, the code will calculate molWeight based on nC, nH and nO
        #2 means [O2 N2]
        self.nSpecies = nSpecies
        self.molWeight = zeros(self.nSpecies)
        #evapFlux out
        self.evapFlux = zeros(self.nSpecies)
        # area and vol
        self.dia = diameter
        self.h = thickness
        self.len = length
        self.vol = pi * self.dia * self.len * self.h
        self.area = pi * (self.dia - 2 * self.h) * self.len

        self.T = T
        self.P = P
        self.Psat = zeros(self.nSpecies)
        self.massDens = zeros(self.nSpecies)
        self.molDens = zeros(self.nSpecies)
        self.massFrac = zeros(self.nSpecies)
        self.molFrac = zeros(self.nSpecies)
        self.volFrac = zeros(self.nSpecies)
        self.concs = zeros(self.nSpecies)
        self.Dvi = zeros(self.nSpecies)
        #air
        self.airMolFrac = zeros(2)
        self.airP = zeros(2)
        self.airMolWeight = array([32.0, 28.0134]) / 1000.
        self.airMassDens = array([1.429, 1.251])

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
        self.vol = pi * self.dia * self.len * self.h
        self.area = pi * (self.dia - 2 * self.h) * self.len
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

    def setAntoine(self, A, B, C):
        """ use Antoine Equation to get saturated vapor pressure
        note the units in Antoine Eqn is mmHg and C
        P = 10^(A-B/(C+T))
        use the system temperature
        """
        A = array(A);B = array(B);C = array(C)
        # transfer from mmHg to Pa
        # A = A0 + log10(101325/760)
        A = A + 2.124903
        # transfer from C to K
        C = C-273.15
        # Antoine's Equation
        self.Psat = 10. ** (A-B / (C + self.T))
        # get air's partial pressure
        self.update()


    def vaporDens(self):
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
        rhovi = self.vaporDens()
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
        drho/dt=A_l/V_l Qi - A_l/V_l (rhoi - sum(Qi)/sum(rhoi))
        y is mass fractional density
        """
        massFracDens = Y[:-1]
        molFracDens = massFracDens / self.molWeight
        self.massFrac = massFracDens / sum(massFracDens)
        self.molFrac = molFracDens / sum(molFracDens)
        ratio = self.area / self.vol
        Qi = self.vaporDiff(Lv=self.dia)
        Q = sum(Qi)

        dhdt = -1. / sum(massFracDens) * Q
        drhodt = -ratio * Qi - ratio * massFracDens * dhdt
        drhodt = append(drhodt, dhdt)
        self.update()
        return drhodt

    def advance(self, t, plotresult=False):
        y0 = self.concs * self.molWeight
        y0 = append(y0, self.h)
        yt = odeint(self.rightSideofODE, y0, t)
        if(plotresult):
            import matplotlib.pyplot as plt
            plt.plot(t, yt)
            plt.show()
        self.h = yt[-1][-1]
        ytt = yt[-1][:-1]
        #        for iii in range(len(ytt)):
        #            if ytt[iii]<0:
        #                ytt[iii]=0.
        molFracDens = ytt / self.molWeight
        self.concs = molFracDens
        self.molFrac = molFracDens / sum(molFracDens)
        self.massFrac = ytt / sum(ytt)

if __name__ == "__main__":
	dia = 0.14 * 10 ** -3
	L = 0.5 * 10 ** -3
	t = 3 * 10 ** -6

	diesel = LiquidFilmCell(nSpecies=7, T=473, diameter=dia, length=L, thickness=t)
	diesel.setCHO(nC=[11, 13, 11, 25, 16, 18, 10],
                  nH=[24, 28, 10, 34, 26, 40, 44],
                  nO=[0, 0, 0, 0, 0, 0, 0])
	print 'diesel components mol weight is', diesel.molWeight #g/mol
	diesel.setMolDens(molDens=[4945.0, 4182.0, 3700.0,
                      3415.0, 2600.0, 2889., 2729.]) #mol/m3
	print 'diesel components mass density is', diesel.massDens #kg/m3
	diesel.setVolFrac(volFrac=[0.05, 0.19, 0.11, 0.25, 0.12, 0.18, 0.10])
	print 'the mol fraction is ', diesel.molFrac
	print 'the mass fraction is ', diesel.massFrac
	print 'the concentrations are ', diesel.concs
	diesel.setAntoine(
                      A=[6.9722, 7.00756, 7.03592, 7.02867, 7.8148, 7.0153, 7.0842],
                      B=[1569.57, 1690.67, 1826.948, 1830.51, 2396.8, 1932.8, 2054],
                      C=[187.7, 174.22, 195.002, 154.45, 199.5736, 137.6, 120.1])
	print 'the saturated vapor pressure is', diesel.Psat
	print 'the O2 and N2 partial pressure is ', diesel.airP[0], diesel.airP[1]
	print 'the O2 and n2 mole fraction are ', diesel.airMolFrac[0], diesel.airMolFrac[1]
	print 'the mol fraction is ', diesel.molFrac
	print 'the mass fraction is ', diesel.massFrac
	print 'the concentrations are ', diesel.concs
	print 'the vapor densities are ', diesel.vaporDens()
	qi = diesel.vaporDiff(Lv=dia)
	#print 'the mass flux out of the interface ',qi
	print 'the initial h is', diesel.h
	print 'start evaporating'
	diesel.advance(arange(0, 0.309, 0.001),True)
	print 'the concentrations are ', diesel.concs
	print 'the vapor densities are ', diesel.vaporDens()
	print 'the new h is', diesel.h
	print '%f percent film left', diesel.h / t









