from evapor import *
import numpy as np
import matplotlib.pyplot as plt

deposit_PS = []
deposit_old = []
deposit_old2 = [] #no phase separation
"""
It seems old2 has very small difference with old
so we don't need turn off phase separation and do
it all over again
"""
for h in [1,3,5,10,15]:
        dia = 1E-3
        L = 5E-3
        initial_film_thickness = h*1E-6
        T = 250
        temperature_in_C = T

        """
        every cycle, if 8Hz, 125ms, however, only 70ms is used for heating
        """
        final_time = 0.07 # s

        diesel = LiquidFilmCell(T=temperature_in_C+273, diameter=dia, length=L, thickness=initial_film_thickness,
                                EVAPORATION=True, CHEMICAL_REACTION=True, PHASE_SEPARATION=True)
        deposit = diesel.deposit

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
        print "Trying DASSL solver"
        diesel.initialize_solver()
        timesteps=linspace(0,final_time,1001)


        #check the residual works (although the initialise_solver above has just done so)
        #import pdb; pdb.set_trace()
        residual = diesel.residual(diesel.t, diesel.y, diesel.dydt)


        diesel_history_array = numpy.zeros((len(timesteps),diesel.nSpecies), numpy.float64)
        diesel_history_concs = numpy.zeros((len(timesteps),diesel.nSpecies), numpy.float64)
        deposit_history_array = numpy.zeros((len(timesteps),diesel.nSpecies), numpy.float64)    
        for step,time in enumerate(timesteps):
            if time>0 : diesel.advance(time)
            # notice these history_arrays now contain MASSES IN MILLIGRAMS!
            diesel_history_array[step] = diesel.y[:diesel.nSpecies] * diesel.SCALE * diesel.molar_masses*1e6
            diesel_history_concs[step] = diesel.concentrations
            deposit_history_array[step] = diesel.y[diesel.nSpecies:] * diesel.SCALE * diesel.molar_masses*1e6
            print "At time t=%g s"%(diesel.t)
            print "Diesel contains (mg)\n", diesel_history_array[step]
            print "Deposit contains (mg)\n", deposit_history_array[step]

        # get deposit using the previous 'all products' method
        depositmass = 0
        for name in 'n-C11(2)  n-C13(3) \
                  Mnphtln(4)  n-C16(5)  C10bnzn(6)  n-C19(7)  n-C21(8)'.split():
            species_index = diesel.species_names.index(name)
            depositmass += diesel_history_array[-1][species_index]
        depositmass = sum(diesel_history_array[-1]) - depositmass
        deposit_old.append(depositmass)
        deposit_PS.append(deposit.total_mass*1e6)

        
        print 'the deposit mass using old method is ',depositmass #mg
        print 'the deposit per area using old method is ', depositmass/diesel.area*1e-6 #mg/mm^2
        print 'final thickness (micrometres)', diesel.thickness * 1e6
        print 'fraction of original thickness', diesel.thickness/diesel.initial_thickness
        print 'the total mass of deposit phase is %g mg' % (deposit.total_mass * 1e6)
        print 'that amount *6500*60 = %g mg' % (deposit.total_mass * 1e6 * 6500 * 60)

# for h in [1,3,5,10,15]:
#         dia = 1E-3
#         L = 5E-3
#         initial_film_thickness = h*1E-6
#         T = 250
#         temperature_in_C = T

#         """
#         every cycle, if 8Hz, 125ms, however, only 70ms is used for heating
#         """
#         final_time = 0.07 # s

#         diesel = LiquidFilmCell(T=temperature_in_C+273, diameter=dia, length=L, thickness=initial_film_thickness,
#                                 EVAPORATION=True, CHEMICAL_REACTION=True, PHASE_SEPARATION=False)
#         deposit = diesel.deposit

#         print 'diesel components molar mass is', diesel.molar_masses # kg/mol
#         print 'diesel components molar density is', diesel.molar_densities # mol/m3
#         print 'diesel components mass density is', diesel.mass_densities # kg/m3
#         print 'the mol fraction is', diesel.mole_fractions
#         print 'the concentrations are', diesel.concentrations
#         print 'the total vapor pressure using K is',sum(diesel.concentrations/diesel.properties.PartitionCoefficientT(diesel.T))*R*diesel.T
#         print 'the saturated vapor pressure is', diesel.Psats
#         print 'the total vapor pressure using Antoine is ',sum(diesel.Psats*diesel.mole_fractions)
#         print 'the total vapor pressure used in the model is ',sum(diesel.vapor_partial_pressures)

#         print 'the vapor densities are ', diesel.get_vapor_mass_densities()
#         qi = diesel.get_evaporative_flux(Lv=dia)
#         print 'the mass flux out of the interface ',qi
#         print 'the initial h is', diesel.thickness
#         print 'initial concentrations are',diesel.concentrations

#         # import pdb; pdb.set_trace(); # pause for debugging
#         print "Trying DASSL solver"
#         diesel.initialize_solver()
#         timesteps=linspace(0,final_time,1001)


#         #check the residual works (although the initialise_solver above has just done so)
#         #import pdb; pdb.set_trace()
#         residual = diesel.residual(diesel.t, diesel.y, diesel.dydt)


#         diesel_history_array = numpy.zeros((len(timesteps),diesel.nSpecies), numpy.float64)
#         diesel_history_concs = numpy.zeros((len(timesteps),diesel.nSpecies), numpy.float64)
#         deposit_history_array = numpy.zeros((len(timesteps),diesel.nSpecies), numpy.float64)    
#         for step,time in enumerate(timesteps):
#             if time>0 : diesel.advance(time)
#             # notice these history_arrays now contain MASSES IN MILLIGRAMS!
#             diesel_history_array[step] = diesel.y[:diesel.nSpecies] * diesel.SCALE * diesel.molar_masses*1e6
#             diesel_history_concs[step] = diesel.concentrations
#             deposit_history_array[step] = diesel.y[diesel.nSpecies:] * diesel.SCALE * diesel.molar_masses*1e6
#             print "At time t=%g s"%(diesel.t)
#             print "Diesel contains (mg)\n", diesel_history_array[step]
#             print "Deposit contains (mg)\n", deposit_history_array[step]

#         # get deposit using the previous 'all products' method
#         depositmass = 0
#         for name in 'n-C11(2)  n-C13(3) \
#                   Mnphtln(4)  n-C16(5)  C10bnzn(6)  n-C19(7)  n-C21(8)'.split():
#             species_index = diesel.species_names.index(name)
#             depositmass += diesel_history_array[-1][species_index]
#         depositmass = sum(diesel_history_array[-1]) - depositmass
#         deposit_old2.append(depositmass)

        
#         print 'the deposit mass using old method is ',depositmass #mg
#         print 'the deposit per area using old method is ', depositmass/diesel.area*1e-6 #mg/mm^2
#         print 'final thickness (micrometres)', diesel.thickness * 1e6
#         print 'fraction of original thickness', diesel.thickness/diesel.initial_thickness
#         print 'the total mass of deposit phase is %g mg' % (deposit.total_mass * 1e6)
#         print 'that amount *6500*60 = %g mg' % (deposit.total_mass * 1e6 * 6500 * 60)

depositold = np.array(deposit_old)
#depositold2 = np.array(deposit_old2)
depositPS = np.array(deposit_PS)


