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


def explodable_range(
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
        return np.logical_and(stars["mass"] >= mmin, stars["mass"] <= mmax)

    return explodability_func


def get_sn_count(
    stars: SSP,
    t0: Quantity["time"],
    t1: Quantity["time"],
    lifetime_func: Callable[[SSP], Quantity["time"]],
    explodability_func: Callable[[SSP], np.bool],
) -> int:
    """
    Get the number of supernovae that have gone off between time t0 and t1.

    :param stars: Stellar Population
    :type stars: SSP
    :param t0: Start time
    :type t0: Quantity["time"]
    :param t1: End time
    :type t1: Quantity["time"]
    :param lifetime_func: Function to calculate stellar lifetimes
    :type lifetime_func: Callable[[SSP], Quantity["time"]]
    :param explodability_func: Function to determine which stars explode as SNe
    :type explodability_func: Callable[[SSP], np.bool]
    :return: Number of supernovae between t0 and t1
    :rtype: int
    """
    lifetimes = lifetime_func(stars)[explodability_func(stars)]
    return np.sum((lifetimes >= t0) & (lifetimes < t1))


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
    Z = stars["metals"]
    logZ = np.log10(np.clip(Z, 7e-5, 3e-2))
    a0 = 10.13 + 0.07547 * logZ - 0.008084 * np.power(logZ, 2)
    a1 = -4.424 - 0.7939 * logZ - 0.1187 * np.power(logZ, 2)
    a2 = 1.262 + 0.3385 * logZ + 0.05417 * np.power(logZ, 2)
    return (
        np.power(
            10,
            a0
            + a1 * np.log10(stars["mass"].to(u.Msun).value)
            + a2 * np.power(np.log10(stars["mass"].to(u.Msun).value), 2),
        )
        * u.yr
    )
