########
## File containing functions relevant for simple planet structure calculations
########

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import SmoothBivariateSpline
from scipy.optimize import fsolve
from scipy.optimize import brentq
import matplotlib.pyplot as plt
from f_constants import *
from microphysics_params import *

###### functions for solid-cores from Fortney et al. 2007

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

def mass_to_radius_solid(mass,Xiron,Xice):

    if (Xiron > 0. ):
        # use iron function
        if (Xice > 0. ):
            print ("error, cannot use Iron and Ice fraction together")

            return -1.

        else:
            rad = iron_function(np.log10(mass),[Xiron,0.])
    else:
        rad = ice_function(np.log10(mass),[Xice,0.])

    return rad


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

#### functions to calculate the I2_I1 integrals, see Owen & Wu (2017) for definition

def integrand1(x,gamma):

    return x*(1./x-1.)**(1./(gamma-1.))

def integrand2(x,gamma):

    return x**2.*(1./x-1.)**(1./(gamma-1.))

def get_I2_I1(DR_Rc,gamma):

    ratio = np.zeros(np.size(DR_Rc))

    for i in range(np.size(DR_Rc)):

        I2 = quad(integrand2,1./(DR_Rc[i]+1.),1.,args=gamma)
        I1 = quad(integrand1,1./(DR_Rc[i]+1.),1.,args=gamma)

        ratio[i] = I2[0]/I1[0]

    return ratio

def get_I2(DR_Rc,gamma):

    I2 = np.zeros(np.size(DR_Rc))

    for i in range(np.size(DR_Rc)):

        I2[i] = quad(integrand2,1./(DR_Rc[i]+1.),1.,args=gamma)[0]

    return I2

def solve_structure(X,Teq,Mcore,Tkh_Myr,Xiron,Xice):

    # this function solves for the structure of the adiabatic envelope of the planet
    # given an envelope mass fraction, i.e. we wish to find the radius of the
    # radiative-convective boundary

    Rcore = mass_to_radius_solid(Mcore,Xiron,Xice) * earth_radius_to_cm

    # use rough analytic formula from Owen & Wu (2017) to guess

    Delta_R_guess = 2.*Rcore * (X/0.027)**(1./1.31) * (Mcore/5.)**(-0.17)

    lg_Delta_Rrcb_guess = np.log10(Delta_R_guess)

    input_args = [X,Teq,Mcore,Tkh_Myr,Xiron,Xice]

    # use log Delta_Rrcb_guess to prevent negative solutions

    lg_D_Rrcb_sol = fsolve(Rrcb_function,lg_Delta_Rrcb_guess,args=input_args)[0]

    Rrcb_sol = 10.**lg_D_Rrcb_sol + Rcore

    # now find f-factor to compute planet radius

    rho_rcb = get_rho_rcb(lg_D_Rrcb_sol,X,Mcore,Teq,Tkh_Myr,Xiron,Xice)

    # now calculate the densities at the photosphere
    Pressure_phot = (2./3. * (G*Mcore*earth_mass_to_g/\
                            (Rrcb_sol**2.*kappa0*Teq**beta)))**(1./(1.+alpha))
    rho_phot_calc = (mu/kb) * Pressure_phot / Teq

    # now find f factor
    H = kb * Teq * Rrcb_sol ** 2. / ( mu * G * Mcore*earth_mass_to_g)

    f = 1. + (H/Rrcb_sol)*np.log(rho_rcb/rho_phot_calc)

    Rplanet = f*Rrcb_sol

    return Rrcb_sol, f, Rplanet, Rrcb_sol-Rcore

def Rrcb_function(lg_D_Rrcb,input_args):

    # we combine equation 4 and 13 from Owen & Wu (2017) to produce a function to solve

    # first evaluate the density at the radiative convective boundary

    # unpack-input arguments
    X=input_args[0]
    Teq = input_args[1]
    Mcore = input_args[2]
    Tkh_Myr = input_args[3]
    Xiron = input_args[4]
    Xice = input_args[5]

    Rcore = mass_to_radius_solid(Mcore,Xiron,Xice) * earth_radius_to_cm

    Rrcb = 10.**lg_D_Rrcb + Rcore

    Delta_R_Rc = 10.**lg_D_Rrcb / Rcore

    rho_core = Mcore * earth_mass_to_g / (4./3.*np.pi*Rcore**3.)

    rho_rcb = get_rho_rcb(lg_D_Rrcb,X,Mcore,Teq,Tkh_Myr,Xiron,Xice)


    I2 = get_I2(np.array([Delta_R_Rc]),gamma)

    cs2 = kb * Teq / mu

    Xguess=3.*(Rrcb/Rcore)**3.*(rho_rcb/rho_core)*(grad_ab * \
              (G * Mcore * earth_mass_to_g)/(cs2 * Rrcb))**(1./(gamma-1.))*I2

    return Xguess - X

def get_rho_rcb(lg_D_Rrcb,X,Mcore,Teq,Tkh_Myr,Xiron,Xice):

    # evaluate the density at the radiative convective boundary - equation
    # 13 from owen & wu

    Rcore = mass_to_radius_solid(Mcore,Xiron,Xice) * earth_radius_to_cm

    Rrcb = 10.**lg_D_Rrcb + Rcore

    Delta_R_Rc = 10.**lg_D_Rrcb / Rcore

    I2_I1 = get_I2_I1(np.array([Delta_R_Rc]),gamma)

    TKh_sec = Tkh_Myr * 1.e6 * year_to_sec

    rho_rcb = (mu / kb) *(I2_I1*64.*np.pi*sigma*Teq**(3.-alpha-beta)\
               *Rrcb*TKh_sec / (3.*kappa0*Mcore*earth_mass_to_g*X))**(1./(1.+alpha))

    return rho_rcb

def evaluate_X(Delta_Rrcb,Teq,Mcore,Tkh_Myr,Xiron,Xice):
    # combining equations 4 and 13 from owen & Wu it is possible to analytically
    # solve for X given the radius of the radiative convective boundary

    Rcore = mass_to_radius_solid(Mcore,Xiron,Xice) * earth_radius_to_cm

    rho_core = Mcore * earth_mass_to_g / (4./3.*np.pi*Rcore**3.)

    cs2 = kb * Teq / mu

    Rrcb = Delta_Rrcb + Rcore

    Delta_R_Rc = Delta_Rrcb / Rcore

    # dimensionless integrals
    I2_I1 = get_I2_I1(np.array([Delta_R_Rc]),gamma)

    I2 = get_I2(np.array([Delta_R_Rc]),gamma)

    TKh_sec = Tkh_Myr * 1.e6 * year_to_sec

    rho_rcb_without_X_term = (mu / kb) *(I2_I1*64.*np.pi*sigma*\
                             Teq**(3.-alpha-beta)*Rrcb*TKh_sec / \
                             (3.*kappa0*Mcore*earth_mass_to_g))**(1./(1.+alpha))

    LHS = 3.*(Rrcb/Rcore)**3.*(rho_rcb_without_X_term/rho_core)*\
    (grad_ab * (G * Mcore * earth_mass_to_g)/(cs2 * Rrcb))**(1./(gamma-1.))*I2

    X = LHS**(1./(1.+1./(1.+alpha)))

    # now find f-factor to compute planet radius, by accounting for the
    # radiative zone on top of the planet

    rho_rcb = get_rho_rcb(np.log10(Delta_Rrcb),X,Mcore,Teq,Tkh_Myr,Xiron,Xice)

    # now calculate the densities at the photosphere
    Pressure_phot = (2./3. * (G*Mcore*earth_mass_to_g/\
                              (Rrcb**2.*kappa0*Teq**beta)))**(1./(1.+alpha))
    
    rho_phot_calc = (mu/kb) * Pressure_phot / Teq

    # now find f factor
    H = kb * Teq * Rrcb ** 2. / ( mu * G * Mcore*earth_mass_to_g)

    f = 1. + (H/Rrcb)*np.log(rho_rcb/rho_phot_calc)

    Rplanet = f*Rrcb


    return X,f,Rplanet

def evaluate_X_rad(Rp,Teq,Mcore,Tkh_Myr,Xiron,Xice):

    # we only use this if the solver returns a adiabatic layer
    # size smaller than a scale height
    # this occurs when the approximation between equations 2 & 3
    # in owen & Wu (2017) breaks down

    # now calculate the densities at the photosphere
    Pressure_phot = (2./3. * (G*Mcore*earth_mass_to_g/(Rp**2.*kappa0*Teq**beta)))**(1./(1.+alpha))
    rho_phot_calc = (mu/kb) * Pressure_phot / Teq

    H = kb * Teq * Rp ** 2. / ( mu * G * Mcore*earth_mass_to_g)

    # now determine mass in radiative layer
    Rcore = mass_to_radius_solid(Mcore,Xiron,Xice) * earth_radius_to_cm

    rho_base = rho_phot_calc * np.exp(Rp/H*(Rp/Rcore-1.))

    Menv = 4.*np.pi*Rcore**2.*H*rho_base

    X = Menv/(Mcore*earth_mass_to_g)

    return X, Rp/Rcore, Rp


def Rp_solver(Rp,Teq,Mcore,Tkh_Myr,Xiron,Xice):

    #this solves for the planet structure to match radii

    input_args=[]
    input_args.append(Rp)
    input_args.append(Teq)
    input_args.append(Mcore)
    input_args.append(Tkh_Myr)
    input_args.append(Xiron)
    input_args.append(Xice)

    Rcore = mass_to_radius_solid(Mcore,Xiron,Xice) * earth_radius_to_cm

    if (Rp < Rcore):
        #error this is not a valid solution for the mass, need to exit
        return -1., -1., -1.

    lg_D_Rrcb_guess = np.log10(Rp-Rcore)
    lg_D_Rrcb_sol = fsolve(Rp_solver_function,lg_D_Rrcb_guess,args=input_args)


    # now evaluate planet structure

    Delta_Rrcb = 10.**lg_D_Rrcb_sol
    #H= kb * Teq * Rp**2. / (mu * G * Mcore * earth_mass_to_g)
    ## Moved this test to elsewhere in the code

    #if (Delta_Rrcb/H < 1.):
        #print("Warning, no convective zone found")
        #X, f, Rplanet = evaluate_X_rad(Rp,Teq,Mcore,Tkh_Myr,Xiron,Xice)
    #else:
    X, f, Rplanet = evaluate_X(Delta_Rrcb,Teq,Mcore,Tkh_Myr,Xiron,Xice)
    
    return X, f, Rplanet


def Rp_solver_function(lg_D_Rrcb,input_args):

    #this function solves for the planet structure to match radii
    Rp=input_args[0]
    Teq = input_args[1]
    Mcore = input_args[2]
    Tkh_Myr = input_args[3]
    Xiron = input_args[4]
    Xice = input_args[5]

    Delta_Rrcb = 10.**lg_D_Rrcb

    # evaluate the envelope mass fraction and planet structure for this guess

    X, f, Rplanet = evaluate_X(Delta_Rrcb,Teq,Mcore,Tkh_Myr,Xiron,Xice)

    return Rp-Rplanet

