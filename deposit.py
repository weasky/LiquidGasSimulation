from __future__ import division,print_function
from numpy import pi,exp, append, arange, array, zeros, linspace, ceil
from evapor import LiquidFilmCell, OneAtm
import matplotlib.pyplot as plt
import copy

class EngineCycleSolver:
    """simulate various engine operating conditions
    """
    def __init__(self,liquidfilm,speed=3000.,stroke=4,run_hours=0.01,sok_time=0.):
        """speed is the engine speed (RPM)
        stroke is the number of strokes per cycle (2 or 4)
        """
        self.speed = speed
        self.run_time = run_hours*3600 # in seconds
        self.secPerCycle = 60.*stroke/speed/2
        self.cycles = ceil(self.run_time / self.secPerCycle)
        self.liquidfilm = liquidfilm
        self.initial_film = copy.deepcopy(liquidfilm)
        self.totalDepositMass = zeros(self.cycles+1)
    
    def getDepositMassPerCycle(self):
        timesteps = linspace(0,self.secPerCycle,101)
        self.liquidfilm.advance(timesteps,plotresult=True)
        liquidMass = self.liquidfilm.concs*self.liquidfilm.molWeight*self.liquidfilm.vol
        depositMass = sum(liquidMass[10:])
        return depositMass

    def getDepositMass(self):
        """the products will participage in future cycle reactions"""
        for cycle in range(self.cycles):
            timesteps = linspace(0,self.secPerCycle,101)
            self.liquidfilm.advance(timesteps,plotresult=True)
            liquidMol = self.liquidfilm.concs*self.liquidfilm.vol
            liquidMass = liquidMol*self.liquidfilm.molWeight
            depositMassInCycle = sum(liquidMass[10:])
            #initial deposit mass is always zero
            #TODO: maybe I should add them together here, using
            #self.totalDepositMass[cycle+1] = self.totalDepositMass[cycle]+depositMassInCycle
            self.totalDepositMass[cycle+1] = depositMassInCycle
            #assume the products' volume is negligible
            self.liquidfilm = copy.deepcopy(self.initial_film)
            self.liquidfilm.concs[10:] = liquidMol[10:]/self.liquidfilm.vol
            if (cycle%50==0 ):
                print('finished cycle %d in total %d cycles.'%(cycle,self.cycles) )
                import time
                time.sleep(2)

        return self.totalDepositMass
            

    def reset(self):
        self.totalDepositMass = zeros(self.cycles+1)
        self.liquidfilm = copy.deepcopy(self.initial_film)


if __name__ == '__main__':
    #params about the nozzle
    dia = 0.14E-3
    L = 0.5E-3
    initial_film_thickness = 3E-6
    run_hours = 100 #run engine 100 hours
    #temperature and pressure
    T=473 #K
    P=OneAtm

    diesel = LiquidFilmCell(T=473,P=1*OneAtm, diameter=dia, length=L, thickness=initial_film_thickness)
    engine = EngineCycleSolver(diesel)
    print ('time per cycle is ',engine.secPerCycle) #3000rpm should give 40ms
    """two deposit models for now:
    1. all lumped products go to deposit layer and they won't join the reaction in the future cycles
    2. all lumped products go to deposit layer and they join the reactions in future cycles
    TODO: I will add deposit 2 layers model and Richard will add phase seperation model"""
    #first model
    depositPerCycle = engine.getDepositMassPerCycle()
    print ('after one cycle, the liquid film has %.2f %%  left.' % (engine.liquidfilm.thickness/initial_film_thickness*100))
    print('the deposit mass formation rate in first cycle is %f g/mm^2.' % (depositPerCycle/engine.liquidfilm.area*1000))
    print('the deposit mass in first cycle is %g g.' % (depositPerCycle*1000 ))
    print('if each cycle is totally independent, the total deposit mass is %f g.' % (depositPerCycle*1000*engine.cycles))
    #second model
    totDepositMass = engine.getDepositMass()
    print(totDepositMass.shape)
    plt.plot(arange(engine.cycles+1),totDepositMass)
    plt.show()
    print('each cycle, the products are persistant, the final deposit mass is %f g.' % (totDepositMass[-1]) )
    
