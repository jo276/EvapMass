########
## File containing functions relevant for mass-loss calculations
########

import numpy as np
import planet_structure as ps
from scipy.optimize import minimize_scalar
from f_constants import *
from scipy.optimize import brentq


# this file contains the details of mass-loss models
def efficiency(Mp,Rp):

    # calculates the mass-loss efficiency

    # constant efficiency

    # eff = 0.1
    # return eff

    # variable efficiency from Owen & Wu (2017)

    vesc = np.sqrt ( 2. * G * Mp / Rp )

    eff = 0.1 * (vesc / 1.5e6)**(-2.)

    return eff

def tmdot_rocky(planet,tmdot_Myr,Xiron,Xice):

    # this function calculates a scaled-mass-loss time for the rocky planet_systems
    # namely it find the envelope mass-fraction at which the mass-loss timescale
    # is maximised

    input_args=[]
    input_args.append(planet.Teq)
    input_args.append(planet.mass)
    input_args.append(tmdot_Myr)
    input_args.append(Xiron)
    input_args.append(Xice)

    DR_min = 0.1*planet.radius * earth_radius_to_cm
    DR_max = 5.*planet.radius * earth_radius_to_cm

    result = minimize_scalar(tmdot_structure,bounds=(DR_min,DR_max),args=input_args,method="bounded")

    if (result.success):

        DR_maximised = result.x
        # now evaluate planet properties
        X,f,Rplanet = ps.evaluate_X(DR_maximised,planet.Teq,planet.mass,tmdot_Myr,Xiron,Xice)

        eff = efficiency(planet.mass*earth_mass_to_g,Rplanet)

        tmdot_max = X * planet.mass **2. * planet.a**2. * eff / Rplanet **3.

        return DR_maximised,X,f,Rplanet,tmdot_max

    else:
        print('Failed to find solution for maximum mass-loss timescale for rocky planet:')
        print(planet.name)
        return -1.


def tmdot_structure(Delta_Rrcb,input_args):

    # this is the mass-loss time-scale function we wish to maximise

    # unpack-input arguments

    Teq = input_args[0]
    Mcore = input_args[1]
    Tkh_Myr = input_args[2]
    Xiron = input_args[3]
    Xice = input_args[4]

    X,f,Rplanet = ps.evaluate_X(Delta_Rrcb,Teq,Mcore,Tkh_Myr,Xiron,Xice)

    # as we wish to maximise the mass-loss timescale at fixed-core mass
    # we wish to find the Delta_Rrcb which maximises:

    # tmdot \propto func_to_max = X * eta (Mp,Rp) / Rp^3

    eta = efficiency(Mcore*earth_mass_to_g,Rplanet)

    func_to_max = X * eta / Rplanet **3.

    # since rountine works to minimize function return 1/func_to_max

    return 1./func_to_max

def find_hardest_rocky(system,tmdot_Myr,Xiron,Xice):

    tmax = 0.
    index_hardest=-1 # set to negative number to spot errors
    Xmax = -1.

    counter = 0
    for planet in system.planets:
        if (planet.rocky_or_gaseous == 1):
            DR_maximised,X,f,Rplanet,t_rocky = tmdot_rocky(planet,tmdot_Myr,Xiron,Xice)
            if (t_rocky > tmax ):
                tmax = t_rocky
                Xmax = X
                index_hardest = counter
        counter += 1

    system.index_rocky_to_scale = index_hardest

    # add mass-loss timescale for rocky planet to hardest one to strip
    system.planets[index_hardest].tmdot = tmax
    system.planets[index_hardest].Xstrip = Xmax


    return system

def min_mass_gaseous(p_rocky,p_gas,Tkh_scale_myr,Xiron,Xice,age_Myr):
    # we wish to find the minimum mass for the gaseous planet given
    # the mass-loss time-scale for the rocky planet

    # the first thing we wish to do is actually check a solution exists
    # we do this by finding the maximum-mass loss timescale and checking it does
    # go above the one we want

    # first find the maximum core mass to check-up to (i.e. when 10% of the
    # planet's radius comes from the convective envelope)

    Rcore = p_gas.radius/1.1

    Mcore_max = ps.solid_radius_to_mass(Rcore,Xiron,Xice)

    input_args=[]
    input_args.append(p_gas.Teq)
    input_args.append(Tkh_scale_myr)
    input_args.append(Xiron)
    input_args.append(Xice)
    input_args.append(p_gas.radius*earth_radius_to_cm)
    input_args.append(age_Myr)
    input_args.append(p_gas.a)
    input_args.append(0.) # tmdot_want is zero for the maximisation

    Mcore_min_try = 0.1 # just use 0.1 earth masses as this is not constrining
    # check solver will actually give a solution for this low a core-mass upto 1 Earth mass

    while Mcore_min_try < Mcore_max:

        sol = tmdot_gas_minimise(Mcore_min_try,input_args)
        if sol < 0.:
            Mcore_min_try += 0.1

            #print("Increased Mcore_min_try to:",Mcore_min_try, "Mcore_max = ", Mcore_max)
        else:
            if (1./sol < p_rocky.tmdot):
                # proceed with this lower bound
                break
            else:
                # solver cannot find a lower core mass that it can solve structure for
                # recommend trying a smaller increase in the min core mass
                # or if happening for 0.1 Mearth then this is a reasonable upper-limit
                #print("Error could not find suitable lower mass bound")

                #return -Mcore_min_try, -5
                return -5., -5

    if sol < 0.:
        # solver cannot find a suitable lower mass bound it can solve for
        #print("Error, could not find a lower bound to solve from")
        return -6., -6

    result = minimize_scalar(tmdot_gas_minimise,bounds=(Mcore_min_try,Mcore_max),args=input_args,method="bounded")

    if result.success:

        tmdot_max_mcore = result.x
        # now evaluate mass-loss timescale at maximum core-mass
        tmdot_max = 1./(tmdot_gas_minimise(tmdot_max_mcore,input_args))

        if (tmdot_max < p_rocky.tmdot):
            # there is no solution
            p_gas.min_mass_sol=False
            p_gas.max_tmdot_core=tmdot_max_mcore

            #return -Mcore_max, -2
            return -2., -2
    else:

        # failed to find maximum core mass for gaseous planet

        return -3., -3

    # now use root-finder between interval to solve for max-mass

    tmdot_min = 1./(tmdot_gas_minimise(Mcore_min_try,input_args))

    if (tmdot_min > p_rocky.tmdot):

        # solution must be lower than Mcore_min_try Mearth

        #return -Mcore_min_try, -1
        return -1., -1

    # to make it this far in the code the minimum mass from photoevaporation
    # must lie between Mcore_min_try Mearth and tmdot_max_mcore so we can use an interval
    # root finder to solve for the value

    # now add actuall
    input_args[7]=(p_rocky.tmdot)

    # work with log mass to prevent negative solutions

    Msol = brentq(mass_gas_to_solve,np.log10(Mcore_min_try),np.log10(tmdot_max_mcore),args=input_args)

    Msol = 10.**Msol

    return Msol, 1

def mass_gas_to_solve(lg_Mcore,input_args):

    Teq = input_args[0]
    Tkh_Myr = input_args[1]
    Xiron = input_args[2]
    Xice = input_args[3]
    Rplanet_now = input_args[4]
    Planet_age_Myr = input_args[5]
    sep_cm = input_args[6]
    tmdot_want = input_args[7]

    # first thing is give the radius of the planet today and current core mass
    # we need to evaluate the envelope mass fraction

    Mcore = 10.**lg_Mcore

    X,f,Rp_solver = ps.Rp_solver(Rplanet_now,Teq,Mcore,Planet_age_Myr,Xiron,Xice)

    # we now want to check that the Rp solver converged by comparing Rp_solver
    # with the actual radius of the planet


    if (np.fabs(Rp_solver-Rplanet_now)/Rplanet_now > 1e-4):
        # Rp_solver failed to converge
        #print("Rp_solver failed to converge, trying again")
        #print(Rp_solver)
        #print(Rplanet_now)
        return -1.

    # now given this envelope mass fraction calculate its radius at time to scale to
    Rrcb_sol, f, Rplanet_t_scale, Delta_R =  ps.solve_structure(X,Teq,Mcore,Tkh_Myr,Xiron,Xice)

    # now check solve_structure converged to correct solution by analytically
    # evaluating the envelope mass fraction given the value of Rrcb

    Xtest,ftest,Rplanet_test = ps.evaluate_X(Delta_R,Teq,Mcore,Tkh_Myr,Xiron,Xice)

    if (np.fabs(Xtest-X)/X > 1e-4):
        # solve structure failed
        #print("solve_structure failed")

        return -1.

    Rplanet = Rplanet_t_scale

    eff = efficiency(Mcore*earth_mass_to_g,Rplanet)

    tmdot_gaseous = X * Mcore**2. * sep_cm**2. * eff / Rplanet**3.

    return tmdot_gaseous - tmdot_want

def tmdot_gas_minimise(Mcore,input_args):

    input_args[7]=0. # makesure tmdot_want is set to zero

    lg_Mcore = np.log10(Mcore)

    # evaluate mass-loss timescale for gaseous planet

    tmdot_gas = mass_gas_to_solve(lg_Mcore,input_args)

    # now return 1 over solution as we want to maximise tmdot_gas

    return 1./tmdot_gas
