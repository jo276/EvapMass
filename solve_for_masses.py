########
## File containing structures and routines to perform min-mass calcualtions
########

import numpy as np
import mass_loss as ml
import planet_structure as ps
from scipy.optimize import fsolve
from f_constants import *

# create objects to hold planet data

class psystem:
    def __init__(self,system_name):
        self.name=system_name
        self.star=star()
        self.number_of_planets=0
        self.planets = []  # list of planets
        self.index_rocky_to_scale = -1 # set to minus one initially for error checking

    def add_planet(self,name,radius,radius_std,period,period_std):
        self.planets.append(planet(name,radius,radius_std,period,period_std))
        self.number_of_planets += 1

    def update_planet_info(self):
        # calculate extra info for planets
        # equilibrium temperature
        for p in self.planets:
            p.a = ( G * self.star.mass * Msun / (2.*np.pi / (p.period * day_to_sec))**2. ) **(1./3.)
            p.Teq = self.star.Teff * np.sqrt(self.star.radius*Rsun/(2.*p.a))

    def above_or_below_valley(self, valley_loc=1.8):
        count = 0
        for p in self.planets:
            p.rocky_or_gaseous = define_valley(p.radius, valley_loc, p.radius_std, p.period)
            count += p.rocky_or_gaseous

        self.num_rocky_planets = count

        return

    def mass_rocky(self,Xiron,Xice):
        # caculate the mass of all the "rocky" planets
        for p in self.planets:
            if (p.rocky_or_gaseous == 1):
                p.mass = ps.solid_radius_to_mass(p.radius,Xiron,Xice)

        return


class star:
    def __init__(self):
        self.mass=0.
        self.mass_std = 0.
        self.age = 5000. # default of 5 Gyr given in Myr
        self.age_std = 0.
        self.radius = 0.
        self.radius_std = 0.
        self.Teff = 0.
        self.Teff_std = 0.

class planet:
    def __init__(self,Name, Radius, Radius_std, Period, Period_std):
        self.name = Name
        self.radius = Radius
        self.radius_std = Radius_std
        self.period = Period
        self.period_std = Period_std
        self.rocky_or_gaseous = -1 # negative to start for error checking
        self.mass = -1.  # set to negative to start for error checking

def setup_systems(planet_systems,Tmdot,Xiron,Xice, valley_loc=1.8):
    # this function takes the list of planetary multi-planet mixed_systems
    # it removes any systems that does not contain both types of planets above
    # and below the gap, it also calculates the mass of the rocky planet from a
    # mass radius-relationship and creates a new list of the mixed systems
    # it also calculates which rocky planet would be hardest to strip.

    # this function takes an input the planet_systems list, but also
    # the mass-loss timescale (age) at which to do this, default is 100Myr (owen & wu 2017)
    # also the core composition can be described in terms of the iron and ice fractions
    # with the fortney et al. 2007 relations only one fraction can be non-zero

    mixed_systems = []
    for system in planet_systems:
        system.update_planet_info()
        system.above_or_below_valley(valley_loc)

        if ((system.num_rocky_planets > 0) and (system.num_rocky_planets < system.number_of_planets) ):

            # system contains both supposed rocky and gaseous planets
            system.mass_rocky(Xiron,Xice)
            # calculate which rocky planet is hardest for evaporation to strip
            system = ml.find_hardest_rocky(system,Tmdot,Xiron,Xice)
            mixed_systems.append(system)

    return mixed_systems

def estimate_min_masses(planet_systems,Tmdot,Xiron,Xice):

    # this function takes the hardest to strip rocky planet and scales all
    # gaseous planets to this one to estimate the mimimum mass required to survive
    # based on the photoevaporation model

    for system in planet_systems:
        if system.index_rocky_to_scale < 0:
            # hardest rocky planet to strip not found:
            print("Error, cannot find which rocky planet to scale from")
            return

        # now loop through all planets in system and evaluate the min mass for
        # those planets flagged as gaseous

        for planet in system.planets:
            if planet.rocky_or_gaseous==0:
                # gaseous planet
                min_mass, converged = ml.min_mass_gaseous(system.planets[system.index_rocky_to_scale],planet,Tmdot,Xiron,Xice,system.star.age)
                planet.min_mass = min_mass
                planet.min_mass_converged = converged # if equal to 1 solution converged, otherwise not

    return

def above_or_below_valley(system, valley_loc):
    for planet in system.planets:
        planet.rocky_or_gaseous = define_valley(planet.radius, valley_loc, \
        planet.radius_std, planet.period)

    return

def define_valley(radius, valley_loc, radius_std, period):
    if radius < valley_loc:
        #rocky planet
        return 1
    else:
        #gaseous planet
        return 0


def error_calculation(value, upper_error, lower_error):
    #function to account for guassian like errors on the stellar and planetary parameters
    list_of_values = []
    value = value
    u = upper_error
    l = lower_error
    for number in range(0,1000,1):
        value_std = np.sqrt(u**2. + l**2)
        value_random = np.random.normal(value, value_std)
        list_of_values.append(value_random)

    new_value = np.mean(list_of_values)
    new_value_std = np.std(list_of_values)

    return new_value , new_value_std
