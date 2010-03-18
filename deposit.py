from __future__ import division,print_function
from numpy import pi,exp, append, arange, array, zeros, linspace, ceil
from evapor import LiquidFilmCell, OneAtm
import matplotlib.pyplot as plt
import copy

class EngineDepositSolver:
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

class TestDepositSolver:
    """simulate various engine operating conditions
    """
    def __init__(self,liquidfilm,air_p,fuel_p,air_pulse,fuel_pulse,rest_pulse,total_pulse,run_hours=20):
        """air_p,fuel_p: pressure in psi
        air_pulse,fuel_pulse,total_pulse: pulse width in ms

        """
        self.air_p = air_p*6894.757 # pascal
        self.fuel_p = fuel_p*6894.757
        self.air_pusle = air_pulse/1000 #seconds
        self.fuel_pulse = fuel_pulse/1000
        self.rest_pulse = rest_pulse / 1000
        self.total_pulse = total_pulse/1000
        self.run_time = run_hours*3600 # in seconds
        self.secPerCycle = (total_pulse -air_pulse - fuel_pulse - rest_pulse)/1000. #heating time
        self.cycles = ceil(self.run_time / self.total_pulse)
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
    initial_film_thickness = 1E-6
    run_hours = 100 #run engine 100 hours
    #temperature and pressure
    T=273+250 #K
    P=OneAtm
    speed = 1500 #RPM,default 3000

    diesel = LiquidFilmCell(T=T,P=P, diameter=dia, length=L, thickness=initial_film_thickness)
    engine = EngineDepositSolver(diesel,speed=speed,run_hours=run_hours)
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
    # #second model
    # totDepositMass = engine.getDepositMass()
    # print(totDepositMass.shape)
    # plt.plot(arange(engine.cycles+1),totDepositMass)
    # plt.show()
    # print('each cycle, the products are persistant, the final deposit mass is %f g.' % (totDepositMass[-1]) )
    """Test configuration simulation"""
    print('======================Start test deposit simulation=========================')
    dia = 1e-3 #1mm
    L = 38e-3 #38 mm in heating, 50mm height tube
    L = 50e-3
    initial_film_thickness = 5e-6 #not sure, may get an estimate from fluid dynamics
    run_hours = 20 #run test 20 hours
    #temperature and pressure
    T=273+250 #K
    P=1.0*OneAtm

    diesel = LiquidFilmCell(T=T,P=P, diameter=dia, length=L, thickness=initial_film_thickness)
    test = TestDepositSolver(liquidfilm=diesel,run_hours=run_hours,
                             air_p=20,fuel_p=25,
                             air_pulse=50,fuel_pulse=2.5,rest_pulse=50,total_pulse=1000)
    #first model
    depositPerCycle = test.getDepositMassPerCycle()
    print ('In test, heating time per cycle is ',test.secPerCycle) 
    print ('after one cycle, the liquid film has %.2f %%  left.' % (test.liquidfilm.thickness/initial_film_thickness*100))
    print('the deposit mass formation rate in first cycle is %f g/mm^2.' % (depositPerCycle/test.liquidfilm.area*1000))
    print('the deposit mass in first cycle is %g g.' % (depositPerCycle*1000 ))
    print('if each cycle is totally independent, the total deposit mass is %f g.' % (depositPerCycle*1000*test.cycles))

    
    
