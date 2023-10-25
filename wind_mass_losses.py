import numpy as np


#########################################################################################
#######################     OB STAR PRESCRIPTIONS MDOT & VWIND     ######################
#########################################################################################

################################
#########  STARBURST99 #########
################################

def mdot_sb99_ostar(Mstar, Age, iso):
    # Equation 1 of Leiherer et al. 1992
    # gives their version of the mass loss rate from
    # main sequence O-stars based on fits to several
    # "simulations" which interatively solve the steady state
    # momentum equation and radiation transfer to get 
    # the radiative force term
    # Mstar : stellar mass in units of Msun
    # Age   : log10(age/year)
    # iso   : MIST Isochrone file to be read from
    # returns log10(Mdot_wind/Msun/year)
 
    # age index for MIST ischrone file
    age_ind = iso.age_index(Age)
    # array of initial masses for isochrone file
    imass = iso.isos[age_ind]['initial_mass']

    # logLbol  : log10(Bolometric luminosity in units of Lsun)
    # value at "right" included to make sure stars that have exploded are not counted
    logLbol = np.interp(Mstar, imass, iso.isos[age_ind]['log_L'], right=-1)
    # logTeff  : log10(Effective Stellar Atmospheric Temperature (K))
    logTeff = np.interp(Mstar, imass, iso.isos[age_ind]['log_Teff'])
    # logZ     : log10(Atmospheric metallicity mass fraction)
    logZ = np.interp(Mstar, imass, iso.isos[age_ind]['log_surf_z'])
    
    # log10(Zsun), taken from MIST Files
    logZsun = -1.8450984742990064
 
    term1 = -24.06 + 2.45*logLbol - 1.1*np.log10(Mstar)
    term2 = 1.31*logTeff + 0.80*(logZ - logZsun)
    result = term1 + term2

    # make sure that exploded stars are not counted
    minval = np.min(result) - 3
    result = (term1 + term2)*(logLbol > 0) + minval*((logLbol < 0))
    return result

def vwind_sb99_ostar(Mstar, Age, iso):
    # Equation 2 of Leiherer et al. 1992
    # gives their version of the wind velocity for
    # main sequence O-stars based on fits to several
    # mehtods and arguments as above
    # returns log10(vwind/km/s)

    age_ind = iso.age_index(Age)
    imass = iso.isos[age_ind]['initial_mass']

    logLbol = np.interp(Mstar, imass, iso.isos[age_ind]['log_L'], right=-1)
    logTeff = np.interp(Mstar, imass, iso.isos[age_ind]['log_Teff'])
    logZ = np.interp(Mstar, imass, iso.isos[age_ind]['log_surf_z'])
    logZsun = -1.8450984742990064

    term1 = 1.23 - 0.3*logLbol + 0.55*np.log10(Mstar)
    term2 = 0.64*logTeff + 0.13*(logZ - logZsun)
    result = term1 + term2

    # make sure that exploded stars are not counted
    minval = np.min(result) - 3
    result = (term1 + term2)*(logLbol > 0) + minval*(logLbol < 0)
    return result

################################
####  VINK ET AL. 2000,2001 ####
################################

def mdot_vink_ostar_useGamma(Mstar, Age, iso):
    # Gives the mass loss prescription adopted in Vink et al. 2001
    # This is laid out most clearly in Section 8 of that paper and
    # depends upon Equations 15 and Equation 23-25
    # this is determined from Monte Carlo radiation transfer solutions
    # Mstar : stellar mass in units of Msun
    # Age   : log10(age/year)
    # iso   : MIST Isochrone file to be read from
    # returns log10(Mdot_wind/Msun/year)

    # age index for MIST ischrone file
    age_ind = iso.age_index(Age)
    # array of initial masses for isochrone file
    imass = iso.isos[age_ind]['initial_mass']

    # logLbol  : log10(Bolometric luminosity in units of Lsun)
    logLbol = np.interp(Mstar, imass, iso.isos[age_ind]['log_L'],right=-1)
    # logTeff  : log10(Effective Stellar Atmospheric Temperature (K))
    logTeff = np.interp(Mstar, imass, iso.isos[age_ind]['log_Teff'])
    # logZ     : log10(Atmospheric metallicity mass fraction)
    logZ = np.interp(Mstar, imass, iso.isos[age_ind]['log_surf_z'])
    
    # log10(Zsun), taken from MIST Files
    logZsun = -1.8450984742990064

    # Effective temperature of transition across the primary 
    # bi-stability jump, given in kiloKelvin (converted to K here)
    # Equations 14 & 15 of Vink et al. 2001
    # this assumes Gamma_e = 0.434 which isn't generally true
    Gamma = 7.66E-5 * 0.325 * (10**logLbol)/Mstar
    charrho = -14.94 + (3.1857 * Gamma) + (0.85 * (logZ - logZsun)) ; 
    TEff_jump = 1e3*(61.2 + 2.59*charrho)

    # Equation 25 of Vink et al 2001 for the mass loss
    # on the cool side of the primary bi-stability jump
    # constant term accounts for vinf/vesc = 1.3 term
    term1_lo = -6.388 + 2.210*(logLbol - 5) - 1.339*np.log10(Mstar/30)
    term2_lo = 1.07*(logTeff - np.log10(2e4)) + 0.85*(logZ-logZsun)

    # Equation 24 of Vink et al 2001 for the mass loss
    # on the hot side of the primary bi-stability jump
    # constant term accounts for vinf/vesc = 2.6 term
    term1_hi = -6.837 + 2.194*(logLbol - 5) - 1.313*np.log10(Mstar/30)
    term2_hi = 0.933*(logTeff - np.log10(4e4)) - 10.92*((logTeff - np.log10(4e4))**2) + 0.85*(logZ-logZsun)

    # Break Results into those meant to be applied above
    # and below the bi-stability jump
    term_lo = (term1_lo + term2_lo)*(10**logTeff<TEff_jump)
    term_hi = (term1_hi + term2_hi)*(10**logTeff>TEff_jump)
    result = term_lo + term_hi

    # make sure that exploded stars are not counted
    minval = np.min(result) - 3
    result = (term_lo + term_hi)*(logLbol > 0) + minval*(logLbol < 0)
    return result

def mdot_vink_ostar(Mstar, Age, iso):
    # Gives the mass loss prescription adopted in Vink et al. 2001
    # This is laid out most clearly in Section 8 of that paper and
    # depends upon Equations 15 and Equation 23-25
    # this is determined from Monte Carlo radiation transfer solutions
    # Mstar : stellar mass in units of Msun
    # Age   : log10(age/year)
    # iso   : MIST Isochrone file to be read from
    # returns log10(Mdot_wind/Msun/year)

    # age index for MIST ischrone file
    age_ind = iso.age_index(Age)
    # array of initial masses for isochrone file
    imass = iso.isos[age_ind]['initial_mass']

    # logLbol  : log10(Bolometric luminosity in units of Lsun)
    logLbol = np.interp(Mstar, imass, iso.isos[age_ind]['log_L'],right=-1)
    # logTeff  : log10(Effective Stellar Atmospheric Temperature (K))
    logTeff = np.interp(Mstar, imass, iso.isos[age_ind]['log_Teff'])
    # logZ     : log10(Atmospheric metallicity mass fraction)
    logZ = np.interp(Mstar, imass, iso.isos[age_ind]['log_surf_z'])
    
    # log10(Zsun), taken from MIST Files
    logZsun = -1.8450984742990064

    # Effective temperature of transition across the primary 
    # bi-stability jump, given in kiloKelvin (converted to K here)
    # Equations 14 & 15 of Vink et al. 2001
    # this assumes Gamma_e = 0.434 which isn't generally true
    TEff_jump = 1e3*(61.2 + 2.59*(-13.636 +0.889*(logZ - logZsun)))

    # Equation 25 of Vink et al 2001 for the mass loss
    # on the cool side of the primary bi-stability jump
    # constant term accounts for vinf/vesc = 1.3 term
    term1_lo = -6.388 + 2.210*(logLbol - 5) - 1.339*np.log10(Mstar/30)
    term2_lo = 1.07*(logTeff - np.log10(2e4)) + 0.85*(logZ-logZsun)

    # Equation 24 of Vink et al 2001 for the mass loss
    # on the hot side of the primary bi-stability jump
    # constant term accounts for vinf/vesc = 2.6 term
    term1_hi = -6.837 + 2.194*(logLbol - 5) - 1.313*np.log10(Mstar/30)
    term2_hi = 0.933*(logTeff - np.log10(4e4)) - 10.92*((logTeff - np.log10(4e4))**2) + 0.85*(logZ-logZsun)

    # Break Results into those meant to be applied above
    # and below the bi-stability jump
    term_lo = (term1_lo + term2_lo)*(10**logTeff<TEff_jump)
    term_hi = (term1_hi + term2_hi)*(10**logTeff>TEff_jump)
    result = term_lo + term_hi

    # make sure that exploded stars are not counted
    minval = np.min(result) - 3
    result = (term_lo + term_hi)*(logLbol > 0) + minval*(logLbol < 0)
    return result


def vwind_vink_ostar(Mstar, Age, iso):
    # Gives the terminal wind velocity of the wind for O stars from
    # Vink et al. 2001 prescirption. It just decides which side of 
    # the bi-stability jump the star is on and then sets 
    # vwind = 1.3, 2.6 * v_esc
    # Mstar : stellar mass in units of Msun
    # Age   : log10(age/year)
    # iso   : MIST Isochrone file to be read from
    # returns log10(vwind/km/s)


    # age index for MIST ischrone file
    age_ind = iso.age_index(Age)
    # array of initial masses for isochrone file
    imass = iso.isos[age_ind]['initial_mass']

    # logTeff  : log10(Effective Stellar Atmospheric Temperature (K))
    logTeff = np.interp(Mstar, imass, iso.isos[age_ind]['log_Teff'],right=-1)
    # Rstar : stellar radius in units of solar masses
    Rstar = 10**np.interp(Mstar, imass, iso.isos[age_ind]['log_R'])
    # logZ     : log10(Atmospheric metallicity mass fraction)
    logZ = np.interp(Mstar, imass, iso.isos[age_ind]['log_surf_z'])
    
    # log10(Zsun), taken from MIST Files
    logZsun = -1.8450984742990064

    # Effective temperature of transition across the primary 
    # bi-stability jump, given in kiloKelvin (converted to K here)
    # Equations 14 & 15 of Vink et al. 2001
    # this assumes Gamma_e = 0.434 which isn't generally true
    TEff_jump = 1e3*(61.2 + 2.59*(-13.636 +0.889*(logZ - logZsun)))

    # escape velocity in km/s
    vesc = 617.7*np.sqrt(Mstar/Rstar)
    # get ratio of wind to escape velocity
    ratio = 1.3*(10**logTeff < TEff_jump) + 2.6*(10**logTeff > TEff_jump)
    result = np.log10(ratio*vesc)

    minvalue = np.min(result) - 3
    result =  np.log10(ratio*vesc)*(logTeff > 0) + minvalue*(logTeff < 0)
    return result

## VINK WEAK WIND

def mdot_vinkww_ostar(Mstar, Age, iso):
    # Gives the mass loss rate for the Vink prescription
    # but tries to account for the "weak wind" problem
    # by taking Mdot -> Mdot/100 if Mstar < 25 Msun
    mdot_vink = mdot_vink_ostar(Mstar, Age, iso)
    return mdot_vink*(Mstar>25) + (mdot_vink-2)*(Mstar<25)

################################
######  VINK & SANDER 2021 #####
################################

def mdot_VS21_ostar(Mstar, Age, iso, use2001=False):
    # Basically a copy of the routines provided as an attachment
    # to the Vink & Sander 2021 paper that was and update to Vink's
    # thesis work. Unlike most implementations of Vink's prescription,
    # this actually takes the Gamma_e dependence into account
    # Mstar : stellar mass in units of Msun
    # Age   : log10(age/year)
    # iso   : MIST Isochrone file to be read from
    # use2001 : boolean to choose a from Z^a relation either 2001 or 2021
    # returns log10(Mdot_wind/Msun/year)

    # Start by reading in all necessary information from MIST 
    # age index for MIST ischrone file
    age_ind = iso.age_index(Age)
    # array of initial masses for isochrone file
    imass = iso.isos[age_ind]['initial_mass']

    # logLbol  : log10(Bolometric luminosity in units of Lsun)
    logLbol = np.interp(Mstar, imass, iso.isos[age_ind]['log_L'],right=-1)
    # logTeff  : log10(Effective Stellar Atmospheric Temperature (K))
    logTeff = np.interp(Mstar, imass, iso.isos[age_ind]['log_Teff'])
    # logZ     : log10(Atmospheric metallicity mass fraction)
    logZ = np.interp(Mstar, imass, iso.isos[age_ind]['log_surf_z'])
    # log10(Zsun), taken from MIST Files
    logZsun = -1.8450984742990064

    # luminosity-to-mass ratio: "Gamma_e"
    Gamma = 7.66E-5 * 0.325 * (10**logLbol)/Mstar

    # Decide which metallicity exponent to use
    dZ = 0.85 if (use2001) else 0.42

    # Calculate Teff of Bi-stability jumps
    # characteristic density for the bi-stability jumps
    # calculated via Eq. (23) from Vink et al. (2001)
    charrho = -14.94 + (3.1857 * Gamma) + (dZ * (logZ - logZsun)) ; 
    # Jump temperatures via Eqs. (15) from Vink et al. (2001) and (6) from Vink et al. (2000)
    T1 = ( 61.2 + (2.59 * charrho) ) * 1000.
    T2 = ( 100. + (6.0 * charrho) ) * 1000.


    Teff = 10**logTeff
    # get ratio of wind-velocity to escape velocity
    ratio = 2.6*(Teff > T1) + 0.7*(Teff < T2) + 1.3*(Teff > T2)*(Teff < T1)
    offset = -6.697*(Teff > T1) + -5.99*(Teff < T2) + -6.688*(Teff > T2)*(Teff < T1)
    
    lograt = np.log10(ratio/2.)
    logT40 = (logTeff - np.log10(4e4))
    logT20 = (logTeff - np.log10(2e4))

    # Below the hotter bi-stability jump, we always use the 2001 exponent
    dZ = 0.85
    logMdot1 = offset + 2.210*(logLbol - 5) - 1.339 * np.log10(Mstar/30.) \
        - 1.601 * lograt + 1.07 * logT20 + dZ * (logZ - logZsun)

    # Above the hotter bi-stability jump, dZ depends on user's choice
    dZ = 0.85 if (use2001) else 0.42
    logMdot2 = offset + 2.194*(logLbol - 5) - 1.313 * np.log10(Mstar/30.) \
        - 1.226 * lograt + 0.933 * logT40 - 10.92 * (logT40**2) + dZ * (logZ - logZsun)

    # apply the correct mass loss in the right regim
    logMdot = logMdot1*(Teff<T1) + logMdot2*(Teff>T1)

    # make sure this is only applied to reasonable regimes and
    # make sure the dead stars are not counted
    select = (T1>T2)*1.0*(logLbol>0)
    minvalue = np.min(logMdot) - 3
    result = logMdot*select + minvalue*(1 - select)

    return result

def vwind_VS21_ostar(Mstar, Age, iso, use2001=False):
    # A very simple version of the wind velocity, using basically what
    # is given in Vink+00,01 but actually taking into account the two
    # bi-stability jumps as in the above mass-loss function
    # Mstar : stellar mass in units of Msun
    # Age   : log10(age/year)
    # iso   : MIST Isochrone file to be read from
    # use2001 : boolean to choose a from Z^a relation either 2001 or 2021
    # returns log10(vwind/km/s)

    # Start by reading in all necessary information from MIST 
    # age index for MIST ischrone file
    age_ind = iso.age_index(Age)
    # array of initial masses for isochrone file
    imass = iso.isos[age_ind]['initial_mass']

    # logLbol  : log10(Bolometric luminosity in units of Lsun)
    logLbol = np.interp(Mstar, imass, iso.isos[age_ind]['log_L'],right=-1)
    # logTeff  : log10(Effective Stellar Atmospheric Temperature (K))
    logTeff = np.interp(Mstar, imass, iso.isos[age_ind]['log_Teff'])
    # logZ     : log10(Atmospheric metallicity mass fraction)
    logZ = np.interp(Mstar, imass, iso.isos[age_ind]['log_surf_z'])
    # Rstar : stellar radius in units of solar masses
    Rstar = 10**np.interp(Mstar, imass, iso.isos[age_ind]['log_R'])
    # log10(Zsun), taken from MIST Files
    logZsun = -1.8450984742990064

    # luminosity-to-mass ratio: "Gamma_e"
    Gamma = 7.66E-5 * 0.325 * (10**logLbol)/Mstar

    # Decide which metallicity exponent to use
    dZ = 0.85 if (use2001) else 0.42

    # Calculate Teff of Bi-stability jumps
    # characteristic density for the bi-stability jumps
    # calculated via Eq. (23) from Vink et al. (2001)
    charrho = -14.94 + (3.1857 * Gamma) + (dZ * (logZ - logZsun)) ; 
    #Jump temperatures via Eqs. (15) from Vink et al. (2001) and (6) from Vink et al. (2000)
    T1 = ( 61.2 + (2.59 * charrho) ) * 1000.
    T2 = ( 100. + (6.0 * charrho) ) * 1000.

    if (np.any(T1 < T2)):
        print("Temperature of Bi-Stability jumps don't make sense...")
        assert(False)

    # escape velocity in km/s
    vesc = 617.7*np.sqrt(Mstar/Rstar)

    Teff = 10**logTeff
    # get ratio of wind-velocity to escape velocity
    ratio = 2.6*(Teff > T1) + 0.7*(Teff < T2) + 1.3*(1.*(Teff > T2) + 1.*(Teff < T1))
    result = np.log10(vesc*ratio)

    # make sure the result is applied in reasonable regimes and
    # make sure dead stars aren't counted
    select = (T1> T2)*1.0*(logLbol>0)
    minvalue = np.min(result) - 3
    result = result*select + minvalue*(1 - select)
    return result


################################
####  BJÖRKLUND ET AL. 2021 ####
################################

def mdot_BJKD_ostar(Mstar, Age, iso):
    # Mass loss rates for OB Stars from the Björklund et al. 2021
    # CMF calculations as quoted from the Appendix of Vink's review
    # interestingly there is no Teff dependence here
    # Mstar : stellar mass in units of Msun
    # Age   : log10(age/year)
    # iso   : MIST Isochrone file to be read from
    # returns log10(Mdot/(Msun/year))

    # Start by reading in all necessary information from MIST 
    # age index for MIST ischrone file
    age_ind = iso.age_index(Age)
    # array of initial masses for isochrone file
    imass = iso.isos[age_ind]['initial_mass']

    # logLbol  : log10(Bolometric luminosity in units of Lsun)
    logLbol = np.interp(Mstar, imass, iso.isos[age_ind]['log_L'],right=-1)
    # logZ     : log10(Atmospheric metallicity mass fraction)
    logZ = np.interp(Mstar, imass, iso.isos[age_ind]['log_surf_z'])
    # log10(Zsun), taken from MIST Files
    logZsun = -1.8450984742990064

    term1 = -5.55 + 0.79*(logZ - logZsun)
    term2 = (logLbol - 6)*(2.16 - 0.32*(logZ - logZsun))
    result = term1 + term2

    # make sure dead stars aren't counted
    minvalue = np.min(result) - 3
    result = result*(logLbol>0) + minvalue*(logLbol<0)
    return result

def vwind_BJKD_ostar(Mstar, Age, iso):
    # wind velocity for OB Stars from the Björklund et al. 2021
    # CMF calculations as quoted from the Appendix of Vink's review
    # derived using the quoted wind momentum relation
    # Mstar : stellar mass in units of Msun
    # Age   : log10(age/year)
    # iso   : MIST Isochrone file to be read from
    # returns log10(Mdot/(Msun/year))

    # Start by reading in all necessary information from MIST 
    # age index for MIST ischrone file
    age_ind = iso.age_index(Age)
    # array of initial masses for isochrone file
    imass = iso.isos[age_ind]['initial_mass']

    # logLbol  : log10(Bolometric luminosity in units of Lsun)
    logLbol = np.interp(Mstar, imass, iso.isos[age_ind]['log_L'],right=-1)
    # logZ     : log10(Atmospheric metallicity mass fraction)
    logZ = np.interp(Mstar, imass, iso.isos[age_ind]['log_surf_z'])
    # logRstar : log10(stellar radius in units of solar masses)
    logRstar = np.interp(Mstar, imass, iso.isos[age_ind]['log_R'])
    # log10(Zsun), taken from MIST Files
    logZsun = -1.8450984742990064

    logMdot = mdot_BJKD_ostar(Mstar, Age, iso)
    term1 = -1*(logMdot + 0.5*logRstar) - 1.55 +0.46*(logZ - logZsun)
    term2 = (logLbol - 6)*(2.07 -0.73*(logZ - logZsun))
    result = term1 + term2

    # make sure dead stars aren't counted
    minvalue = np.min(result) - 3
    result = result*(logLbol>0) + minvalue*(logLbol<0)
    return result

#########################################################################################
#####################     WOLF-RAYET PRESCRIPTIONS MDOT & VWIND     #####################
#########################################################################################

def mdot_NL00_WR(Mstar, Age, iso):
    # Mass loss rates for Wolf-Rayet stars as given by 
    # Nugis & Lamers 2000 based on empirical calibrations
    # Mstar : stellar mass in units of Msun
    # Age   : log10(age/year)
    # iso   : MIST Isochrone file to be read from
    # returns log10(Mdot/(Msun/year))
 
    # age index for MIST ischrone file
    age_ind = iso.age_index(Age)
    # array of initial masses for isochrone file
    imass = iso.isos[age_ind]['initial_mass']

    # logLbol  : log10(Bolometric luminosity in units of Lsun)
    logLbol = np.interp(Mstar, imass, iso.isos[age_ind]['log_L'],right=-1)
    # Y3    : Atmospheric He^3 mass fraction
    Y3 = np.interp(Mstar, imass, iso.isos[age_ind]['surface_he3'])
    # Y4    : Atmospheric He^3 mass fraction
    Y4 = np.interp(Mstar, imass, iso.isos[age_ind]['surface_he4'])
    # logY  : log10(surface mass fraction of Helium)
    logY = np.log10(Y3+Y4)
    # logZ     : log10(Atmospheric metallicity mass fraction)
    logZ = np.interp(Mstar, imass, iso.isos[age_ind]['log_surf_z'])

    # Zsun in their assumed model
    Zsun = 0.018
    term1 = -11 + 1.29*logLbol
    term2 = 1.73*logY + 0.47*(logZ - np.log10(Zsun))
    result = term1 + term2

    # make sure that dead stars are not counted
    minvalue = np.min(result) - 3
    result = (term1 + term2)*(logLbol>0) + minvalue*(logLbol<0)
    return result

def mdot_SV20_WR(Mstar, Age, iso):
    # Mass loss rate of Sander & Vink 2020 taken from the 
    # appendix of Vink's 2021 Review. Based on a large 
    # numer of theoretical CMF calculations using the POWR code
    # Mstar : stellar mass in units of Msun
    # Age   : log10(age/year)
    # iso   : MIST Isochrone file to be read from
    # returns log10(Mdot/(Msun/year))

    # age index for MIST ischrone file
    age_ind = iso.age_index(Age)
    # array of initial masses for isochrone file
    imass = iso.isos[age_ind]['initial_mass']

    # logLbol  : log10(Bolometric luminosity in units of Lsun)
    logLbol = np.interp(Mstar, imass, iso.isos[age_ind]['log_L'],right=-1)
    # this mass loss prescription is really driven by the 
    # iron abundance and so shouldn't account for the self-enriched
    # metallicity but simply the initial iron abundance relative to solar
    logZ = iso.abun['[Fe/H]']

    # luminosity-to-mass ratio: "Gamma_e"
    Gamma = 7.66E-5 * 0.325 * (10**logLbol)/Mstar
    
    Gammaeb = -0.324*logZ + 0.244
    c = -0.44*logZ + 9.15
    logMdot_off = 0.23*logZ - 2.61

    result = 2.932*np.log10(-1*np.log10(1-Gamma)) - np.log10(2)*((Gammaeb/Gamma)**c) + logMdot_off

    # make sure that dead stars are not counted
    minvalue = np.min(result) - 3
    result = (result)*(logLbol>0) + minvalue*(logLbol<0)
    return result

def vwind_SLUG_WR(Mstar, Age, iso, mdot_func = mdot_NL00_WR):
    # Following the SLUG implementation we assume that the 
    # WR momentum input is given by Lbol/c and use a given Mdot
    # function to then infer the wind velocity
    # Mstar : stellar mass in units of Msun
    # Age   : log10(age/year)
    # iso   : MIST Isochrone file to be read from
    # mdot_func : function to get mass loss rates taking (Mstar, Age, iso)
    # returns : log10(v_wind in km/s)
 
    # age index for MIST ischrone file
    age_ind = iso.age_index(Age)
    # array of initial masses for isochrone file
    imass = iso.isos[age_ind]['initial_mass']

    # logLbol  : log10(Bolometric luminosity in units of Lsun)
    logLbol = np.interp(Mstar, imass, iso.isos[age_ind]['log_L'],right=-1)

    # Lbol/c in units of Msun km/s/year
    L_over_c = 2.036e-8*(10**logLbol)
    mdot = 10**mdot_func(Mstar, Age, iso)
    result = np.log10(L_over_c/mdot)

    # make sure that dead stars are not counted
    minvalue = np.min(result) - 3
    result = (np.log10(L_over_c/mdot))*(logLbol>0) + minvalue*(logLbol<0)
    return result