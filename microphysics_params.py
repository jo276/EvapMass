###### Thermodynamic microphysics parameters

from f_constants import *


## opacity from form kappa = kappa0*P^alpha*T^beta
alpha = 0.68 # pressure dependence of opacity
beta = 0.45 # temperature dependence of opacity
kappa0 = 10**(-7.32) # opacity constant


mu = 2.35 * mh # solar metallicity gas
gamma = 5./3. # ratio of specific rho_earth_cgs

grad_ab = (gamma-1.)/gamma
