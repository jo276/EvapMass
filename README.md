
EvapMass uses the photoevaporation driven evolution model of Owen & Wu (2013,2017) to predict the minimum masses of mini-Neptunes in multi-planet systems to be consistent with the photoevaporation model. If planetary masses are known it can also be used to test the system for photoevaporation - if the minimum masses predicted are consisted with the measured masses of the mini-Neptunes, then the system is consistent with photoevaporation.

The planetary system requires both a planet above and below the radius valley to be useful for this test. Where the valley is defined is dependent on the system's host star - make a careful and informed choice. We recommend checking the radius valley slopes predicted by van Eylen et al. (2018,2021).

There is an example jupyter notebook for the Kepler-36 system (Intro_notebook.ipynb) and a notebook which can be easily used to test any system (general_notebook.ipynb) (Feb 2023 update).

The paper on this code can be read here: https://academic.oup.com/mnras/article/491/4/5287/5663631


Update - August 2022

A. Error fix - the efficiency parameter was incorrectly implemented in the mass_loss.py routines for the escape velocity scaling. Thanks goes to Ryan Cloutier and Madison Van Wyngarden for noticing this. This does not change the results shown in Owen & Campos Estrada (2020).

B. The mass-loss efficiency parameter can now be chosen with the 'eff_option' flag:

Option 1 - constant efficiency
Option 2 - escape velocity scaling
Option 3 - Use of the mass-loss rates from Owen & Jackson (2012)

C. New error flag - If the code detects that the planet's atmosphere (for the solved minimum) mass detects that the size of the convection zone is smaller than one-scale height it changes the error_flag to 2. Care should be taken when considering these planet structures as the approximation between Equations 2 and 3 in Owen & Wu (2017) breaks down for this case.

Update - February 2023

The location of the valley can now be chosen by the user with the flag 'valley_loc'. This is in Earth radii. The default is set to 1.8 Rearth, which is generally only true for FGK stars. So be careful if you are testing an M-dwarf system.
