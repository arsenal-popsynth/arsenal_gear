==========================
The Arsenal Parameter File
==========================


The arsenal parameter file is organized into blocks related to different physical modules
in the code. We describe each block below. Some of these choices may be more or less
compatabile with one another. Where there is an egregious example of poorly matched
choices, the code will provide a warning.

Population Parmeters
********************

This block specifies properties of the population of stars that you would like to
simulate. This includes the distribution of stars as a function of mass (the IMF) as well
as distributions over properties of binaries within the population. Right now the code
only allows for a single metallicity and age for the poulation or a "simple stellar
population" (SSP).

Choices related to these options are detailed below.

Stellar Evolution
*****************

This block controls options for the stellar evolution code you would like to use for both
single and binary star evolution. These choices essentially map the inital properties of
a star (i.e. its mass for a single star or binary properties for a binary) to things like
its mass, luminosity, and effective temperature as a funciton of its age.

Choices related to these parameters are detailed below.

Atmosphere Modeling
*******************

This block determines how the properties of a star (its luminosity, mass, radius,
effective temperature, etc.) lead to radiative emission from its atmosphere.

Mechanical Feedback
*******************

This block gives choices for how stars return mass and energy through both stellar winds
(of multiple types) and supernovae (of multiple types). Included within these blocks are
choices for the yields of each feedback mechanism or the mass fraction of each element
that is returned.

