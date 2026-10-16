"""
be_data_structures.py
==========

This file defines the data structures used to store results from binary
evolution models. They build on those in stellar_evolution/se_data_structures.py:
each star is a StellarTrack, and single stars from binary evolution codes use
StellarTrack and TrackSet directly.
"""

from dataclasses import dataclass

import astropy.constants as const
import astropy.units as u
import numpy as np
from astropy.units import Quantity

from ..stellar_evolution.se_data_structures import StellarTrack, TrackSet


def _kepler_semimajor_axis(
    m_tot: Quantity["mass"], period: Quantity["time"]
) -> Quantity["length"]:
    """Semi-major axis of an orbit with total mass m_tot and the given period."""
    a3 = const.G * m_tot * period**2 / (4 * np.pi**2)
    return np.cbrt(a3).to(u.au)


@dataclass
class BinaryTrack:
    """
    Data class for storing one binary system: two StellarTracks sharing an age
    grid, plus the system's initial orbital parameters. The evolution of the
    orbit (period and eccentricity over time) is not stored yet, see #45.

    Attributes:
        primary: track of the initially more massive star (primary.m0 = M1)
        secondary: track of its companion (secondary.m0 = M2 = q * M1)
        q0: initial mass ratio q = M2 / M1
        P0: initial orbital period
        e0: initial orbital eccentricity
    """

    primary: StellarTrack
    secondary: StellarTrack
    q0: float
    P0: Quantity["time"]
    e0: float

    def __post_init__(self):
        if len(self.primary.age) != len(self.secondary.age):
            raise ValueError("primary and secondary tracks must share an age grid")

    @property
    def age(self) -> Quantity["time"]:
        """Ages of the binary track points, shared by both stars."""
        return self.primary.age

    @property
    def mass1(self) -> Quantity["mass"]:
        """Current masses of the primary star."""
        return self.primary.mass

    @property
    def mass2(self) -> Quantity["mass"]:
        """Current masses of the secondary star."""
        return self.secondary.mass

    @property
    def teff1(self) -> Quantity["temperature"]:
        """Effective temperatures of the primary star."""
        return self.primary.teff

    @property
    def teff2(self) -> Quantity["temperature"]:
        """Effective temperatures of the secondary star."""
        return self.secondary.teff

    @property
    def lbol1(self) -> Quantity["power"]:
        """Bolometric luminosities of the primary star."""
        return self.primary.lbol

    @property
    def lbol2(self) -> Quantity["power"]:
        """Bolometric luminosities of the secondary star."""
        return self.secondary.lbol

    @property
    def rad1(self) -> Quantity["length"]:
        """Radii of the primary star."""
        return self.primary.rad

    @property
    def rad2(self) -> Quantity["length"]:
        """Radii of the secondary star."""
        return self.secondary.rad

    @property
    def grav1(self) -> Quantity["acceleration"]:
        """Surface gravities of the primary star."""
        return self.primary.grav

    @property
    def grav2(self) -> Quantity["acceleration"]:
        """Surface gravities of the secondary star."""
        return self.secondary.grav

    @property
    def mass_ratio(self) -> np.ndarray[np.float64]:
        """Current mass ratios M2 / M1 of the binary track points."""
        return (self.mass2 / self.mass1).decompose().value

    @property
    def a0(self) -> Quantity["length"]:
        """Initial orbital semi-major axis, from Kepler's third law."""
        return _kepler_semimajor_axis(self.primary.m0 + self.secondary.m0, self.P0)


@dataclass
class BinaryTrackSet(TrackSet):
    """
    Data class for storing a set of binary evolution tracks. Inherited
    attributes describe the systems: masses are the initial primary masses,
    min_ages/max_ages each system's age range, and the variable names apply
    to both the primary and the secondary.

    Attributes (in addition to TrackSet's):
        tracks: list of BinaryTrack objects
        q0: initial mass ratio of each system
        P0: initial orbital period of each system
        e0: initial orbital eccentricity of each system
    """

    tracks: list[BinaryTrack]
    q0: np.ndarray[np.float64]
    P0: Quantity["time"]
    e0: np.ndarray[np.float64]

    @property
    def a0(self) -> Quantity["length"]:
        """Initial orbital semi-major axis of each system, from Kepler's third law."""
        return _kepler_semimajor_axis(self.masses * (1 + self.q0), self.P0)
