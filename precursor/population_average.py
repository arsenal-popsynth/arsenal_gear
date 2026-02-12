# Functions for doing Population Synthesis

import numpy as np
from artpop.stars import sample_imf


def gen_mass_samp(Mstar_tot, Nlist):
    # function to get a sample of stars from a Kroupa IMF
    # that have about a total mass of Mstar_tot, in units
    # of solar masses
    # output is list of Nlist # of list of masses in solar units

    # average mass for the Kroupa IMF, calculated separately
    mavg = 0.57113728

    # list of mass lists to create
    mass_samples = []
    for i in range(Nlist):
        Nsamp = int(1.2 * Mstar_tot / mavg)
        mtrial = sample_imf(
            Nsamp,
            m_min=0.01,
            m_max=200,
            imf={"a": [-0.3, -1.3, -2.3], "b": [0.08, 1.0]},
        )
        while np.sum(mtrial) < Mstar_tot:
            mtrial = sample_imf(
                Nsamp,
                m_min=0.01,
                m_max=200,
                imf={"a": [-0.3, -1.3, -2.3], "b": [0.08, 1.0]},
            )
        try:
            nsel = np.where(np.cumsum(mtrial) < Mstar_tot)[0][-1]
            toss = np.random.uniform()
            # toss a coin to decide to go slightly over or slightly under target
            if toss > 0.5:
                mass_samples.append(mtrial[: nsel + 1])
            else:
                mass_samples.append(mtrial[: nsel + 2])
        except:
            print(np.where(np.cumsum(mtrial) < Mstar_tot)[0], np.sum(mtrial), Mstar_tot)

    return mass_samples


def gen_Lbol(ms, iso_mist, age):
    # gives the summed bolometric luminosity for a mass sample of stars "ms"
    # prescription based on interpolating masses to "iso_mist" isochrone
    # ms: list of masses in solar masses
    # iso_mist: MIST isochrone file
    # age: log10(age/year) has to be between 5.0 and 1.3 in 0.05 increment
    # returns in units of solar luminosities

    # get siochrone information from MIST File
    age_ind = iso_mist.age_index(age)
    logLbol = iso_mist.isos[age_ind]["log_L"]
    imass = iso_mist.isos[age_ind]["initial_mass"]

    # interpolate onto the masses being considered
    lLbs = np.interp(ms, imass, logLbol)
    return np.sum(10**lLbs)


def get_Lbol_evol(ms_list, iso_mist, ages):
    # for a given set of mass_samples (ms_list) this generates the
    # Bolometric luminosity evolution over the list of ages for a given
    # isochrone (iso_mist) and returns the median and
    # interquartile range of the evolution
    # ms_list     : list of mass arrays in units of Msun
    # iso_mist    : theoretical MIST ischrone file to interpolate from
    # ages        : list of ages in log10(age/year)
    # returns in units of solar luminosities

    # loop over ages and call gen_Lbol() function for each
    (Lbol, Lbol_25, Lbol_75) = ([], [], [])
    for age in ages:
        # get total array of total Mdot for each mass sample
        Lbolsum_list = np.array([gen_Lbol(ms, iso_mist, age) for ms in ms_list])
        # derive median and interquartile range over independent mass samples
        Lbol.append(np.median(Lbolsum_list))
        Lbol_25.append(np.quantile(Lbolsum_list, 0.25))
        Lbol_75.append(np.quantile(Lbolsum_list, 0.75))

    # return in array format
    return (np.array(Lbol_25), np.array(Lbol), np.array(Lbol_75))


def gen_mdot(ms, iso_mist, age, mdot_func, op="OB"):
    # gives the sum of mass loss rates according to the "mdot_func"
    # prescription based on interpolating masses to "iso_mist" isochrone
    # ms: list of masses in solar masses
    # iso_mist: MIST isochrone file
    # age: log10(age/year) has to be between 5.0 and 1.3 in 0.05 increment
    # mdot_func: function to provide log10(Mdot) takes (Mstar, Age, iso)
    # op : string, either "OB" or "WR" selects which stars should be passed
    #      onward for mass loss calculation

    # get siochrone information from MIST File
    age_ind = iso_mist.age_index(age)
    logTeff = iso_mist.isos[age_ind]["log_Teff"]
    imass = iso_mist.isos[age_ind]["initial_mass"]
    surfX = iso_mist.isos[age_ind]["surface_h1"]

    # interpolate onto the masses being considered
    lTs = np.interp(ms, imass, logTeff)
    sXs = np.interp(ms, imass, surfX)

    if op == "OB":
        # definition of what an Ostar is from Leitherer et al. 1992
        # log(Teff)> 3.9 and surface_X > 0.4
        # also select stars with M > 8 Msun
        ssel = np.intersect1d(np.where(lTs > 3.9), np.where(sXs > 0.4))
        ssel = np.intersect1d(np.where(ms > 8)[0], ssel)
    elif op == "WR":
        # use the Leitherer et al. 1992 definition since MIST
        # definition apparently has overlap with the MS
        # also select stars with M > 8 Msun
        ssel = np.intersect1d(np.where(lTs > 4.4), np.where(sXs < 0.4))
        ssel = np.intersect1d(np.where(ms > 8)[0], ssel)
    else:
        print("Called with star option that isn't supported...")
        assert False

    if len(ssel) == 0:
        return 0.0
    else:
        return np.sum(10 ** mdot_func(ms[ssel], age, iso_mist))


def get_mdot_evol(ms_list, iso_mist, ages, mdot_func, op="OB"):
    # for a given set of mass_samples (ms_list) this generates the
    # mass loss evolution history over the list of ages for a given
    # mass loss precription (mdot_func) and returns the median and
    # interquartile range of the evolution
    # ms_list     : list of mass arrays in units of Msun
    # iso_mist    : theoretical MIST ischrone file to interpolate from
    # ages        : list of ages in log10(age/year)
    # mdot_func : function to provide log10(Mdot) for OB Stars takes (Mstar, Age, iso)
    # op : string, either "OB" or "WR" selects which stars should be passed
    #      onward for mass loss calculation

    # loop over ages and call gen_mdot() function for each
    (Mdot, Mdot_25, Mdot_75) = ([], [], [])
    for age in ages:
        # get total array of total Mdot for each mass sample
        Mdotsum_list = np.array(
            [gen_mdot(ms, iso_mist, age, mdot_func, op) for ms in ms_list]
        )
        # derive median and interquartile range over independent mass samples
        Mdot.append(np.median(Mdotsum_list))
        Mdot_25.append(np.quantile(Mdotsum_list, 0.25))
        Mdot_75.append(np.quantile(Mdotsum_list, 0.75))

    # return in array format
    return (np.array(Mdot_25), np.array(Mdot), np.array(Mdot_75))


def gen_Lwind(ms, iso_mist, age, mdot_func, vwind_func, op="OB"):
    # gives the sum of mechanical wind luminosites according to
    # the "mdot_func" and "vwind_func"
    # prescription based on interpolating masses to "iso_mist" isochrone
    # ms: list of masses in solar masses
    # iso_mist: MIST isochrone file
    # age: log10(age/year) has to be between 5.0 and 1.3 in 0.05 increment
    # mdot_func: function to provide log10(Mdot) takes (Mstar, Age, iso)
    # vwind_func : function to provide v_wind for OB Stars takes (Mstar, Age, iso)
    # op : string, either "OB" or "WR" selects which stars should be passed
    #      onward for mass loss calculation
    # returns in units of Lwind/(Msun/year (km/s)^2)

    # get siochrone information from MIST File
    age_ind = iso_mist.age_index(age)
    logTeff = iso_mist.isos[age_ind]["log_Teff"]
    imass = iso_mist.isos[age_ind]["initial_mass"]
    surfX = iso_mist.isos[age_ind]["surface_h1"]

    # interpolate onto the masses being considered
    lTs = np.interp(ms, imass, logTeff)
    sXs = np.interp(ms, imass, surfX)

    if op == "OB":
        # definition of what an Ostar is from Leitherer et al. 1992
        # log(Teff)> 3.9 and surface_X > 0.4
        # also select stars with M > 8 Msun
        ssel = np.intersect1d(np.where(lTs > 3.9), np.where(sXs > 0.4))
        ssel = np.intersect1d(np.where(ms > 8)[0], ssel)
    elif op == "WR":
        # use the Leitherer et al. 1992 definition since MIST
        # definition apparently has overlap with the MS
        # also select stars with M > 8 Msun
        ssel = np.intersect1d(np.where(lTs > 4.4), np.where(sXs < 0.4))
        ssel = np.intersect1d(np.where(ms > 8)[0], ssel)
    else:
        print("Called with star option that isn't supported...")
        assert False

    if len(ssel) == 0:
        return 0.0
    else:
        mdots = 10 ** mdot_func(ms[ssel], age, iso_mist)
        vwinds = 10 ** vwind_func(ms[ssel], age, iso_mist)
        return np.sum(0.5 * mdots * (vwinds**2))


def get_Lwind_evol(ms_list, iso_mist, ages, mdot_func, vwind_func, op="OB"):
    # for a given set of mass_samples (ms_list) this generates the
    # wind luminosity evolution over the list of ages for a given
    # mass loss precription (mdot_func) and wind prescription (vwind_func)
    # and returns the median and interquartile range of the evolution
    # returns in units of solar masses per year * (km/s)^2
    # ms_list     : list of mass arrays in units of Msun
    # iso_mist    : theoretical MIST ischrone file to interpolate from
    # ages        : list of ages in log10(age/year)
    # mdot_func : function to provide log10(Mdot) for OB Stars takes (Mstar, Age, iso)
    # vwind_func : function to provide v_wind for OB Stars takes (Mstar, Age, iso)
    # op : string, either "OB" or "WR" selects which stars should be passed
    #      onward for mass loss calculation

    # loop over ages and call gen_Lwind() function for each
    (Lwind, Lwind_25, Lwind_75) = ([], [], [])
    for age in ages:
        # get total array of total Lwind for each mass sample
        Lwindsum_list = np.array(
            [gen_Lwind(ms, iso_mist, age, mdot_func, vwind_func, op) for ms in ms_list]
        )
        # derive median and interquartile range over independent mass samples
        Lwind.append(np.median(Lwindsum_list))
        Lwind_25.append(np.quantile(Lwindsum_list, 0.25))
        Lwind_75.append(np.quantile(Lwindsum_list, 0.75))

    # return in array format
    return (np.array(Lwind_25), np.array(Lwind), np.array(Lwind_75))


def gen_NSN(ms, iso_mist, age):
    # returns the number of SNe that have exploded for a set of stars with
    # intitial masses "ms" at "age" for isochrone "iso_mist"
    # ms       : mass array in units of Msun
    # iso_mist : theoretical MIST ischrone file to interpolate from
    # age     : age in log10(age/year)

    # get siochrone information from MIST File
    age_ind = iso_mist.age_index(age)
    imass = iso_mist.isos[age_ind]["initial_mass"]

    return 1.0 * len(np.intersect1d(np.where(ms > imass[-1]), np.where(ms > 8)))


def get_NSN_evol(ms_list, iso_mist, ages):
    # for a given set of mass_samples (ms_list) this counts the number of
    # SNe that have occured as a function of time (in "ages" list) and
    # returns the median and interquartile range of the evolution
    # ms_list     : list of mass arrays in units of Msun
    # iso_mist    : theoretical MIST ischrone file to interpolate from
    # ages        : list of ages in log10(age/year)

    # loop over ages and call gen_mdot() function for each
    (NSN, NSN_25, NSN_75) = ([], [], [])
    for age in ages:
        # get total array of total Mdot for each mass sample
        Mdotsum_list = np.array([gen_NSN(ms, iso_mist, age) for ms in ms_list])
        # derive median and interquartile range over independent mass samples
        NSN.append(np.mean(Mdotsum_list))
        NSN_25.append(np.quantile(Mdotsum_list, 0.25))
        NSN_75.append(np.quantile(Mdotsum_list, 0.75))

    # return in array format
    return (np.array(NSN_25), np.array(NSN), np.array(NSN_75))


def gen_MSN(ms_list, iso_mist, m_remnant=1.4):
    # function to attempt to get the mass lost from each SNe
    # basically looks for the last mass before the star exploded
    # and subtracts m_remnant from that
    # ms_list     : list of mass arrays in units of Msun
    # iso_mist    : theoretical MIST ischrone file to interpolate from
    # m_remnant   : assumed remnant mass for all stars in solar masses

    ages = [round(x, 2) for x in np.array(iso_mist.ages)]

    # arrays to return that give the explosion masses and ages for each SN
    (Mw_list, Mexp_list) = ([], [])
    # loop through mass samples
    for ms in ms_list:
        # get sorted massive stars that will supernova
        # in decreasing mass order
        Ms = np.sort(ms[np.where(ms > 8)])[::-1]
        nms = len(Ms)
        # pointer to the star under consideration (we want to march
        # downward in mass to avoid looping through multilpe times)
        mp = 0
        # arrays to add the SN explosion masses and
        # estimate ages at which they explode
        Mws = np.zeros(np.array(ages).shape)
        Mexps = np.zeros(np.array(ages).shape)
        for k, age in enumerate(ages):
            # age index for MIST ischrone file
            age_ind = iso_mist.age_index(age)
            # array of initial masses for isochrone file
            imass = iso_mist.isos[age_ind]["initial_mass"]

            # enter loop of determining SN mass loss if the star under
            # consideration is more massive than the maximum initial mass
            while (mp < nms) and (Ms[mp] > imass[-1]):
                # load properties of isochrones at last age available
                old_age_ind = iso_mist.age_index(ages[k - 1])
                imass_old = iso_mist.isos[old_age_ind]["initial_mass"]
                star_mass = iso_mist.isos[old_age_ind]["star_mass"]
                # EEP = iso_mist.isos[old_age_ind]['EEP']

                # interpolate final mass of star and subtract remnant mass
                # to get mass added back to the ISM
                mstar_end = np.interp(Ms[mp], imass_old, star_mass)
                Mexps[k] += mstar_end - m_remnant
                Mws[k] += Ms[mp] - mstar_end
                mp += 1
                # enter another loop in order to avoid re-loading
                # the "old "isochrone data
                while (mp < nms) and (Ms[mp] > imass[-1]):
                    # but do the same thing
                    mstar_end = np.interp(Ms[mp], imass_old, star_mass)
                    Mexps[k] += mstar_end - m_remnant
                    mp += 1
        Mexp_list.append(np.cumsum(Mexps))
        Mw_list.append(np.cumsum(Mws))
    # an array containing arrays for each cumulative mass evolution
    # as a function of time for each mass list in ms_list
    Mexp_list = np.array(Mexp_list)
    Mw_list = np.array(Mw_list)

    # returns median cumulative mass evolution in time
    return (np.array(ages), np.median(Mexp_list, axis=0), np.median(Mw_list, axis=0))
