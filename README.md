The uses the photoevaporation driven evolution model of Owen & Wu (2013,2017) to predict the minimum masses of planets in multi-planet systems to be consistent with the photoevaporation model. 

The planetary system requires both a planet above and below the radius gap to be useful for this test. 

There is an example jupyter notebook for the Kepler-36 system. 

The paper on this code can be read here: https://academic.oup.com/mnras/article/491/4/5287/5663631


Update - August 2022

A. Error fix - the efficiency parameter was incorrectly implemented in the mass_loss.py routines for the escape velocity scaling. Thanks goes to Ryan Cloutier and Madison Van Wyngarden for noticing this. 

B. The efficiency parameter options have been extended through the use of an eff_option flag:

Option 1 - constant efficiency
Option 2 - escape velocity scaling
Option 3 - Use of the mass-loss rates from Owen & Jackson (2012)

C. New error flag - If the code detects that the planet's atmosphere (for the solved minimum) mass detects that the size of the convection zone is smaller than one-scale height it changes the error_flag to 2. Care should be taken when considering these planet structures as the approximation between Equations 2 and 3 in Owen & Wu (2017) breaks down for this case. 
