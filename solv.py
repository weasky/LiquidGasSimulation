##
# BASED ON: ctml_writer.py
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
# Richard West 2009 edited by Amrit

import math, pylab
import quantities as pq


# import all the classes which write XML
# and all the constant, functions, etc.
#reload(ctml_writer)
from ctml_writer import *
import ctml_writer
# now modify the ones we want to do more

class outsideValidRange(Exception):
	"""was not within valid range for expression"""
	
class species(ctml_writer.species):
	"""a species"""
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

class NASA(ctml_writer.NASA):
	"""NASA polynomial representation of thermo"""
	def ValidTemperature(self,temperature):
		if temperature<self._t[0]: return False
		if temperature>self._t[1]: return False
		return True

	def getThermo(self,T):
		""" getThermo(T) returns (HeatCapacityOverR, EnthalpyOverRT, EntropyOverR) for a given temperature T
		    Raises outsideValidRange exception if T is not within range of polynomial"""
		if not self.ValidTemperature(T): raise outsideValidRange
		if self._pref > 0.0: raise Exception("not sure what to do with customised standard state pressure")
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
		EntropyOverR = c[0]*math.log(T) + c[1]*T + c[2]*T*T/2 + c[3]*T*T*T/3 + c[4]*T*T*T*T/4 + c[6]
		return(HeatCapacityOverR,EnthalpyOverRT,EntropyOverR)
		
	def getGibbsFreeEnergy(self,T):
		"""getGibbsFreeEnergy(T) returns  FreeEnergyOverR at temperature T"""
		(HeatCapacityOverR,EnthalpyOverRT,EntropyOverR)=self.getThermo(T)
		FreeEnergyOverRT=EnthalpyOverRT-EntropyOverR
		return FreeEnergyOverRT*T
		

pq.UnitQuantity('kilocalories', pq.cal*1e3, symbol='kcal')
pq.UnitQuantity('kilojoules', pq.J*1e3, symbol='kJ')
pq.UnitQuantity('kilomoles', pq.mol*1e3, symbol='kmol')
_uConc=pq.Quantity(1,ctml_writer._umol)/ pq.Quantity(1,ctml_writer._ulen)**3
_uTime=pq.Quantity(1,ctml_writer._utime)



class reaction(ctml_writer.reaction):
	"""a chemical reaction"""
	
	def getForwardRateCoefficient(self,T):
		"""returns the forward rate coefficient at a given temperature assuming 
		   kf = [A, n, E] and kf=A (T/1K)^n exp(-E/RT)"""
		kf=self._kf
		
		forwardRateCoefficient=kf[0] / _uConc**(self.getReactantNu()-1) / _uTime
		if kf[1]:
			forwardRateCoefficient *= T**kf[1]
		assert ctml_writer._ue=='kcal/mol', 'unit convrsion not yet implemented'
		R=0.0019872065 # kcal/mol/K
		R=pq.constants.R 
		forwardRateCoefficient*=math.exp(-kf[2]/(R*T))
		print "%s forward k = %s"%(self._e,forwardRateCoefficient)
		return forwardRateCoefficient
		
	def getReverseRateCoefficient(self,T):
		"""returns the reverse rate coefficient at a given temperature"""
		if not self.rev: return 0
 # reaction not reversible
		forwardRateCoefficient=self.getForwardRateCoefficient(T)
		reverseRateCoefficient=forwardRateCoefficient/self.getEquilibriumConstant(T)
		print "%s reverse k = %s"%(self._e,reverseRateCoefficient)
		return reverseRateCoefficient
		
	def getForwardRate(self,T,concentrations):
		"""returns the forward rate of progress, given a temperature and a dictionary of species concentrations"""
		forwardRate=self.getForwardRateCoefficient(T)
		for speciesName,order in self._rxnorder.items():
		#	print "forwardRate*=concentrations[%s]**order = %s**%g"%(speciesName,concentrations[speciesName],order)
			forwardRate=forwardRate*concentrations[speciesName]**order
		print "%s forward rate = %s"%(self._e,forwardRate.simplified)
		return forwardRate
	
	def getReverseRate(self,T,concentrations):
		"""returns the reverse rate of progress, given a temperature and a dictionary of species concentrations"""
		if not self.rev: return 0*_uConc/pq.s # reaction not reversible
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
		"""returns the change in gibbs free energy *over R* for the rxn at a given T"""
		deltaGOverR=0
		for speciesName,order in self._p.items(): # add products
			deltaGOverR += order* getSpeciesByName(speciesName).getGibbsFreeEnergy(T)
		for speciesName,order in self._r.items(): # subtract reactants
			deltaGOverR -= order* getSpeciesByName(speciesName).getGibbsFreeEnergy(T)		
		return deltaGOverR
		
	def getEquilibriumConstant(self,T):
		"""returns the equilibrium constant at a given temperature
		    Keq = exp(-DeltaG/RT) """
		
		#if not sum(self._p.values())==sum(self._rxnorder.values()): print "equilibrium constant calculation currently assumes forward and reverse reactions have same reaction order"
		
		dGoverR=self.getDeltaG(T) # returns delta G over R
		Keq= math.exp(-dGoverR/T) 
		
		Kc=Keq * _uConc ** self.getDeltaNu()
		if self.getDeltaNu(): print "this is wrong. Kc should not equal %s"%Kc
		# the units are right, but it's not _uConc that you should be multiplying by
		
		return Kc
		
		
	
	def getReactantNu(self):
		"""returns the stoichiometry in the forwards direction
		i.e. the number of reactant molecules
		uses self._reactantNu to cache the answer"""
		if hasattr(self,'_reactantNu'): return self._reactantNu
		reactantNu=0
		for speciesName,order in self._r.items(): # add reactants
			reactantNu += order
		self._reactantNu=reactantNu
		return reactantNu
	def getProductNu(self):
		"""returns the stoichiometry in the reverse direction. 
		i.e. the number of product molecules
		uses self._productNu to cache the answer"""
		if hasattr(self,'_productNu'): return self._productNu
		productNu=0
		for speciesName,order in self._p.items(): # add products
			productNu += order
		self._productNu=productNu
		return productNu
		
	def getDeltaNu(self):
		"""returns the change in stoichiometry of the reaction, delta Nu.
		   deltaNu= productNu - reactantNu
		   uses self._deltaNu to cache the answer"""
		if hasattr(self,'_deltaNu'): return self._deltaNu
		self._deltaNu=self.getProductNu()-self.getReactantNu()
		return self._deltaNu


_species=ctml_writer._species
_reactions=ctml_writer._reactions
_speciesnames=ctml_writer._speciesnames

def getSpeciesByName(name):
	for s in ctml_writer._species:
		if s._name==name:
			return s
	return None

def getNetRatesOfCreation(T,concs):
	nrocs=dict.fromkeys(_speciesnames)
	for speciesName in _speciesnames:
		nrocs[speciesName]=pq.Quantity(0.0,'mol/m**3/s')
	for r in _reactions:
		rate=r.getNetRate(T,concs)
		for speciesName,order in r._p.items(): # create products
			nrocs[speciesName]+=rate.simplified*order
			print "rate of creation of %s += %s. Now nroc=%s"%(speciesName,rate.simplified*order,nrocs[speciesName])
		for speciesName,order in r._r.items(): # consume reactants
			nrocs[speciesName]-=rate.simplified*order
			print "rate of consumption of %s += %s. Now nroc=%s"%(speciesName,rate.simplified*order,nrocs[speciesName])
	return nrocs

class FuelComponent():
	""" shouldn't really be part of the solv package. specific to the fuel model."""
	def __str__(self):
		return "<Species %s>"%(self.name)
	def __init__(self,name="",initialVolFraction=0,composition=dict(C=0,H=0,O=0),Antoine=dict(A=0,B=0,C=0),liquidMolarDensity=3500 ):
		self.name=name
		self.initialVolFraction=initialVolFraction
		self.composition=composition
		self.Antoine=Antoine
		self.liquidMolarDensity=liquidMolarDensity * pq.mol / pq.m**3  # mol/m3
		self.initialConcentration=self.liquidMolarDensity*initialVolFraction
		
		
if __name__ == "__main__":
	#reload(ctml_writer)
	import sys, os
	if len(sys.argv)>1:
		file = sys.argv[1]
	else:
		print "using default file, chem.cti"
		file = 'chem.cti' # default input file
	base = os.path.basename(file)
	root, ext = os.path.splitext(base)
	dataset(root)
	if not _species:
		execfile(file)
	
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
	
	T=430 # kelvin
	
	timestep=1e-6 * pq.s #seconds
	
	concsHistory=[[concs[spn].copy() for spn in _speciesnames]]
	
	
	for t in range(5):
		nrocs=getNetRatesOfCreation(T,concs)
		for speciesName,nroc in nrocs.items(): 
			concs[speciesName]+=nroc*timestep  # Simple Euler
			
		concsHistory.append([concs[spn].copy() for spn in _speciesnames])
		
	charray=pylab.array(concsHistory)
	pylab.semilogy(charray)
	pylab.show()
	
		
			


