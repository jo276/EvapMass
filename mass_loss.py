import numpy as np

# constants
G = 6.67408e-8
earth_mass_to_g = 5.972e27
earth_radius_to_cm = 6.371e8

# this file contains the details of mass-loss models
def efficiency(planet):

    # calculates the mass-loss efficiency

    # constant efficiency

    # eff = 0.1
    # return eff

    # variable efficiency from Owen & Wu (2017)

    vesc = np.sqrt ( 2. * G * planet.mass * earth_mass_to_g
    / planet.radius / earth_radius_to_cm )

    eff = 0.1 * (vesc / 1.5e6)**(-2.)

    return eff

def tmdot_rocky(planet):

    # this function calculates a scaled-mass-loss time for the rocky planet_systems
    # this rocky planet with the largest mass-loss timescale is the hardest to
    # remove a H/He atmosphere from and the one we scale gaseous planets to.
    # here we care about the maximum mass-loss efficiency when X = X2 or Delta_R=1
    # in the notation of Owen & Wu (2017)

    #power-law indexes - replace with calculations from microphysics parameters
    n_p = 1.41
    n_M = 1.42
    n_eff = -1.

    eff = efficiency(planet)

    tmdot_scale = eff**(n_eff) * (planet.period)**n_p * (planet.mass)**n_M

    return tmdot_scale

def find_hardest_rocky(system):

    tmax = 0.
    index_hardest=-1 # set to negative number to spot errors

    counter = 0
    for planet in system.planets:
        if (planet.rocky_or_gaseous == 1):
            t_rocky = tmdot_rocky(planet)
            if (t_rocky > tmax ):
                tmax = t_rocky
                index_hardest = counter
        counter += 1

    system.index_rocky_to_scale = index_hardest

    return system
