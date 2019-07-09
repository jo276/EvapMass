import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import solve_for_masses as em

data = pd.read_csv('cks_physical_merged.csv')
planet_systems = []

system_counter = 0
new_system = True

for index, row in data.iterrows():
    if (row['koi_disposition'] != 'FALSE POSITIVE'):
        if (index == 0):
            planet_systems.append(em.psystem(row['id_kic']))
        else:
            if (planet_systems[system_counter].name != row['id_kic']):
                if (planet_systems[system_counter].number_of_planets == 1):
                    # remove planet from list and don't update system_counter
                    del planet_systems[system_counter]
                else:
                    system_counter +=1
                new_system = True
                planet_systems.append(em.psystem(row['id_kic']))

        if (new_system):
            planet_systems[system_counter].star.mass = row['iso_smass']
            planet_systems[system_counter].star.radius = row['iso_srad']
            planet_systems[system_counter].star.Teff = row['iso_steff']
            planet_systems[system_counter].star.age = row['iso_sage'] \
            * 1e9 / 1e6 # want age in Myr
            new_system = False


        # add planet
        planet_systems[system_counter].add_planet(row['id_koicand'],\
        row['iso_prad'],row['koi_period'])

Xiron = 1./3.
Xice = 0.
Tmdot_Myr = 100.

mixed_systems = em.setup_systems(planet_systems,Tmdot_Myr,Xiron,Xice)
em.estimate_min_masses(mixed_systems,Tmdot_Myr,Xiron,Xice)

counter = 0
total_planets=0
for system in mixed_systems:
    print("System is", system.planets[0].name)
    for planet in system.planets:
        if planet.rocky_or_gaseous == 0:
            print(planet.min_mass_converged)

            if (planet.min_mass_converged == -5):
                    counter +=1


            total_planets +=1
            print(planet.min_mass)
            print(planet.radius)

print('Number of inconsistent planets')
print(counter)
print('Fraction is')
print(counter*1./(total_planets*1.))
