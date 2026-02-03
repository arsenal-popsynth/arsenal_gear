"""
sn
=========

Functions for supernovae
"""

from typing import Callable

import astropy.units as u
import numpy as np
from astropy.units import Quantity

from ..population import SSP


def constant_energy(stars: SSP, energy: Quantity["energy"]) -> Quantity["energy"]:
    """
    Returns a constant supernova explosion energy for all stars.

    :param stars: Stellar Population
    :type stars: SSP
    :param energy: Constant explosion energy to return
    :type energy: Quantity["energy"]
    :return: Explosion energy for each star in the population
    :rtype: Quantity["energy"]
    """
    return np.full(len(stars), energy)


def explodable_mass_range(
    mmin: Quantity["mass"], mmax: Quantity["mass"]
) -> Callable[[SSP], np.bool]:
    """
    Returns a function that determines which stars are explodable based on a simple mass range.

    :param mmin: Minimum mass for explodability
    :type mmin: Quantity["mass"]
    :param mmax: Maximum mass for explodability
    :type mmax: Quantity["mass"]
    :return: Function that takes a SSP and returns a boolean array indicating explodability
    """

    def explodability_func(stars: SSP) -> np.bool:
        return np.logical_and(stars >= mmin, stars < mmax)

    return explodability_func


def explodable_lifetime_range(
    tmin: Quantity["time"],
    tmax: Quantity["time"],
    lifetime_func: Callable[[SSP], Quantity["time"]],
) -> Callable[[SSP], np.bool]:
    """
    Returns a function that determines which stars are explodable based on a simple lifetime range.

    :param tmin: Minimum lifetime for explodability
    :type tmin: Quantity["time"]
    :param tmax: Maximum lifetime for explodability
    :type tmax: Quantity["time"]
    :param lifetime_func: Function to calculate stellar lifetimes
    :type lifetime_func: Callable[[SSP], Quantity["time"]]
    :return: Function that takes a SSP and returns a boolean array indicating explodability
    """

    def explodability_func(stars: SSP) -> np.bool:
        lifetimes = lifetime_func(stars)
        return np.logical_and(lifetimes >= tmin, lifetimes < tmax)

    return explodability_func


def lifetimes_Raiteri(stars: SSP) -> Quantity["time"]:
    """
    Stellar lifetimes calculated from
    `Raiteri+ 1996 <https://ui.adsabs.harvard.edu/abs/1996A%26A...315..105R/abstract>`__
    Equation 3.

    :param stars: Stellar Population
    :type stars: SSP
    :return: lifetime of each star in the population
    :rtype: Quantity["time"]
    """
    Z = stars.metals
    logZ = np.log10(np.clip(Z, 7e-5, 3e-2))
    a0 = 10.13 + 0.07547 * logZ - 0.008084 * np.power(logZ, 2)
    a1 = -4.424 - 0.7939 * logZ - 0.1187 * np.power(logZ, 2)
    a2 = 1.262 + 0.3385 * logZ + 0.05417 * np.power(logZ, 2)
    return (
        np.power(
            10,
            a0
            + a1 * np.log10(stars.to(u.Msun).value)
            + a2 * np.power(np.log10(stars.to(u.Msun).value), 2),
        )
        * u.yr
    )


def massloss_Raiteri(stars: SSP) -> Quantity["mass"]:
    """
    Mass loss calculated from
    `Raiteri+ 1996 <https://ui.adsabs.harvard.edu/abs/1996A%26A...315..105R/abstract>`__
    Equation 6.

    :param stars: Stellar Population
    :type stars: SSP
    :return: mass lost by each star in the population
    :rtype: Quantity["mass"]
    """
    return 0.7682 * u.Msun * (stars / u.Msun) ** 1.056


def metals_Raiteri(stars: SSP) -> dict:
    """
    Iron mass loss calculated from
    `Raiteri+ 1996 <https://ui.adsabs.harvard.edu/abs/1996A%26A...315..105R/abstract>`__
    Equation 7 and 8, with total metallicity derived from Asplund+ 2009
    <https://ui.adsabs.harvard.edu/abs/2009ARA%26A..47..481A/abstract> abundances.

    :param stars: Stellar Population
    :type stars: SSP
    :return: dictionary of mass losses by element
    :rtype: dict
    """
    iron = 2.802e-4 * u.Msun * (stars / u.Msun) ** 1.864
    oxygen = 4.586e-4 * u.Msun * (stars / u.Msun) ** 2.721
    metals = 1.06 * iron + 2.09 * oxygen
    return {"Fe": iron, "O": oxygen, "metals": metals}
