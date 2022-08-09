import numpy as np

#from EvapMass.mass_loss import efficiency

get_shape = np.load("eff_shape_file.npy", allow_pickle=True).item()

def get_scaled_hydro_eff(Rp,Mp):

    ## returns the scaled efficiency factor of the Owen & Jackson (2012) evaporation rates

    # Rp and Mp should be in cm and grams

    # the shape is scaled via the escape velocity to a mass of 7.603262769401823e+27 g

    Mp_scale = 7.603262769401823e+27

    ## Radius scale first 

    scaled_eff = 10.**get_shape(np.log10(Rp*Mp_scale/Mp))

    ## this scaled efficiency returns either a scaled value (scaled to max value in table) for the
    ## efficiency or extropolates the efficiency at a constant value of large planets
    ## that would be undergoing Roche Lobe overflow in the Owen & Jackson (2012) tables
    ## E.g. grey region of Figure 5 in Owen & Jackson (2012)


    ### scale mass

    mass_scale = (1. + (np.sqrt(Mp/1e29))**10.)**(1./10.)

    return scaled_eff * mass_scale



