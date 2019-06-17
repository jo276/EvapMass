import numpy as np
import mass_loss as ml
from scipy.optimize import fsolve

# constants
Tsun = 5778. #K
Rsun = 6.9551e10 # cm
Msun = 1.9885e33
au_to_cm = 1.496e13
Period_Earth = 365.25 #days
day_to_sec = 86400.
Teq_Earth = Tsun * np.sqrt(Rsun/(2.*au_to_cm))
G = 6.67408e-8

# create objects to hold planet data

class psystem:
    def __init__(self,system_name):
        self.name=system_name
        self.star=star()
        self.number_of_planets=0
        self.planets = []  # list of planets

    def add_planet(self,name,radius,period):
        self.planets.append(planet(name,radius,period))
        self.number_of_planets += 1

    def update_planet_info(self):
        # calculate extra info for planets
        # equilibrium temperature
        for p in self.planets:
            p.a = ( G * Msun / (2.*np.pi / (p.period * day_to_sec))**2. ) **(1./3.)
            p.Teq = self.star.Teff * np.sqrt(self.star.radius*Rsun/(2.*p.a))

    def above_or_below_valley(self):
        count = 0
        for p in self.planets:
            p.rocky_or_gaseous = define_valley(p.radius,p.period)
            count += p.rocky_or_gaseous

        self.num_rocky_planets = count

        return

    def mass_rocky(self,Xiron,Xice):
        # caculate the mass of all the "rocky" planets
        for p in self.planets:
            if (p.rocky_or_gaseous == 1):
                p.mass = solid_radius_to_mass(p.radius,Xiron,Xice)

        return


class star:
    def __init__(self):
        self.mass=0.
        self.age = 5.e9 # default of 5 Gyr
        self.radius = 0.
        self.Teff = 0.

class planet:
    def __init__(self,Name,Radius,Period):
        self.name = Name
        self.radius = Radius
        self.period = Period
        self.rocky_or_gaseous = 0
        self.mass = -1. # set to negative to start

def setup_systems(planet_systems,Xiron,Xice):
    # this function takes the list of planetary multi-planet mixed_systems
    # it removes any systems that does not contain both types of planets above
    # and below the gap, it also calculates the mass of the rocky planet from a
    # mass radius-relationship and creates a new list of the mixed systems
    # it also calculates which rocky planet would be hardest to strip.

    mixed_systems = []

    for system in planet_systems:
        system.update_planet_info()
        system.above_or_below_valley()

        if ((system.num_rocky_planets > 0) and (system.num_rocky_planets < system.number_of_planets) ):
            # system contains both supposed rocky and gaseous planets
            system.mass_rocky(Xiron,Xice)
            # calculate which rocky planet is hardest for evaporation to strip
            system = ml.find_hardest_rocky(system)
            mixed_systems.append(system)

    return mixed_systems



def above_or_below_valley(system):
    for planet in system.planets:
        planet.rocky_or_gaseous = define_valley(planet.radius,planet.period)

    return

def define_valley(radius,period):
    # use crude single value
    if radius < 1.8:
        return 1
    else:
        return 0

def solid_radius_to_mass(radius,Xiron,Xice):
    # convert a radius into mass using mass-radius relationship

    # use Fortney et al. 2007, mass-radius relationship

    if (Xiron > 0. ):
        # use iron function
        if (Xice > 0. ):
            print ("error, cannot use Iron and Ice fraction together")

            return -1.

        else:
            mass = fsolve(iron_function,np.log10(5.),args=[Xiron,radius])
    else:
        mass = fsolve(ice_function,np.log10(5.),args=[Xice,radius])

    return 10.**mass


def ice_function(lg_mass,inputs):

    X=inputs[0]
    radius = inputs[1]

    R = (0.0912*X + 0.1603) * (lg_mass)**2.
    R += (0.3330*X + 0.7387) * (lg_mass)
    R += (0.4639*X + 1.1193)

    return R-radius

def iron_function(lg_mass,inputs):

    Xiron=inputs[0]
    radius = inputs[1]

    X = 1. - Xiron # rock mass fraction

    R = (0.0592*X + 0.0975) * (lg_mass)**2.
    R += (0.2337*X + 0.4938) * (lg_mass)
    R += (0.3102*X + 0.7932)

    return R-radius
