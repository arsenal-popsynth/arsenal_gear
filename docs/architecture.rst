==============================
Arsenal Gear Code Architecture
==============================


Top-Level Population Class
**************************

The main way that users should interact with the code in python scripts is through the
top level ``StellarPopulation`` class. It has member classes for both single and binary
stellar sub-populations::
    StellarPopulation
    params
        samp_mode              (stochastic, integrative)
        mstar                  (total stellar mass)
        metallicity            (stellar metallicity)
        fbin                   (total binary fraction)
    member classes
        SingleStarPop
        SingleStellarEvolution
        BinaryPop
        BinaryStellarEvolution
        Atmosphere
        MechanicalFeedback
    methods
        wind_mass_loss(times)
        radiation(times, band_data)
        sn_rate(times)
        yields(times)
        save_outputs(times)

Abstract Stellar Population Classes
***********************************

These classes contain information on the initial population of single stars and binaries
as well as ::

    SingleStarPop
    params
        IMF               (kroupa, chabrier, etc.)

    BinaryPop
    params
        IMF               (IMF of the primary)
        fbin_dist         (binary fraction distribution)
        q_dist            (mass ratio distribution)
        P_dist            (period distribution)
        e_dist            (eccentricity distribution)

Abstract Stellar Evolution Classes
**********************************

These classes contain interfaces to stellar evolution codes::

    SingleStellarEvolution
    params
        SingleStarPop
    methods
        props(times)          (get stellar properties)
    
    BinaryStellarEvolution
    params
        BinaryPop
    methods
        props_i(times)
        e(times)              (eccentricity evolution)
        P(times)              (Period evolution)

The  ``props`` methods will return ``StellarProperties`` structures which are containers
for all of the properies of a star that are important for its radiation and feedback.

Abstract Atmosphere Classes
***************************

This provides the radiative output of a star based on its properties such as bolometric
luminosity, mass, radius, effective temperature, metallicity::

    Atmosphere
    params
        StellarProperties
    methods
        get_band(times, band_data)
        get_spec(times)

Abstract Mechanical Feedback Classes
************************************

Classes related to the return of mass, energy, and elements to the ISM::

    MechanicalFeedback
    params
        inputs
    member classes
        Winds
        Supernovae
        Yields

    Winds
    params
        StellarProperties
    methods
        mdotw(times)        (wind mass loss rates)
        vwind(times)        (wind velocities)
    
    Supernovae
    params
        StellarEvolution    (not necessary depending on choices)
    methods
        get_mmax(times)     (get maximum mass still alive)

    Yields
    params
        wind_ops
        sn_ops
    methods
        yi_wind             (wind yields)
        yi_sn(m)            (supernova yields as a funciton of mass)
