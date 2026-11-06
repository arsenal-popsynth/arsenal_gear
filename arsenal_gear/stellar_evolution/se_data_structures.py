"""
se_data_structures.py
==========

This file defines various data structures used in the interpolation of
results from stellar evolution models
"""

from dataclasses import dataclass

import astropy.units as u
import numpy as np
from astropy.units import Quantity


@dataclass
class Isochrone:
    """
    Data class for storing isochrone data.

    Attributes:
        age: Age of the isochrone
        eep_name: variable name for equivalent evolutionary phase,
                  None for models without EEPs
        m0_name: variable name for initial stellar masses
        lteff_name: variable name for log_10(effective temperature/K)
        llbol_name: variable name for log_10(bolometric luminosities/L_sun)
        lrad_name: variable name for log_10(stellar radii/R_sun)
        lgrav_name: variable name for log_10(surface gravities/cm/s^2)
        qs: Dictionary of isochrone quantities. These should minimally include
            EEP, initial_mass, log_Teff, and log_L but can include others.
    """

    age: Quantity["time"]
    eep_name: str | None
    m0_name: str
    lteff_name: str
    llbol_name: str
    lrad_name: str
    lgrav_name: str
    # list of surface abundances, if available
    elems: list[str]
    qs: dict

    @property
    def eep(self) -> np.ndarray[int]:
        """Equivalent evolutionary phase numbers of the isochrone points."""
        if self.eep_name is None:
            raise ValueError("isochrone has no EEPs")
        return self.qs[self.eep_name].astype(int)

    @property
    def m0(self) -> Quantity["mass"]:
        """Initial masses of the isochrone points."""
        return self.qs[self.m0_name] * u.Msun

    @property
    def teff(self) -> Quantity["temperature"]:
        """Effective temperatures of the isochrone points."""
        return np.power(10, self.qs[self.lteff_name]) * u.K

    @property
    def lbol(self) -> Quantity["power"]:
        """Bolometric luminosities of the isochrone points."""
        return np.power(10, self.qs[self.llbol_name]) * u.Lsun

    @property
    def rad(self) -> Quantity["length"]:
        """Radii of the isochrone points."""
        return np.power(10, self.qs[self.lrad_name]) * u.Rsun

    @property
    def grav(self) -> Quantity["acceleration"]:
        """Surface gravities of the isochrone points."""
        return np.power(10, self.qs[self.lgrav_name]) * u.cm / u.s**2


# TODO(@ltancas): think of the best way to make age units consistent
#                 between the isochrone and stellar track data structures
@dataclass
class StellarTrack:
    """
    Data class for storing stellar track data.

    Attributes:
        m0: Initial mass of the stellar track
        eeps: Row numbers (1-based) of the primary EEPs along the track,
              None for models without EEPs
        age_name: variable name for age (in years by default)
        mass_name: variable name for current stellar mass/M_sun,
                   None if not available
        lteff_name: variable name for log_10(effective temperature/K)
        llbol_name: variable name for log_10(bolometric luminosities/L_sun)
        lrad_name: variable name for log_10(stellar radii/R_sun)
        lgrav_name: variable name for log_10(surface gravities/cm/s^2)
        qs: Dictionary of stellar track quantities. These should minimally include
            age, log_Teff, and log_L but can include others.
    """

    m0: Quantity["mass"]
    eeps: np.ndarray[int] | None
    age_name: str
    mass_name: str | None
    lteff_name: str
    llbol_name: str
    lrad_name: str
    lgrav_name: str
    # list of surface abundances, if available
    elems: list[str]
    qs: dict

    @property
    def age(self) -> Quantity["time"]:
        """Ages of the stellar track points."""
        return self.qs[self.age_name] * u.yr

    @property
    def mass(self) -> Quantity["mass"]:
        """Current masses of the stellar track points."""
        if self.mass_name is None:
            raise ValueError("stellar track has no current-mass column")
        return self.qs[self.mass_name] * u.Msun

    @property
    def teff(self) -> Quantity["temperature"]:
        """Effective temperatures of the stellar track points."""
        return np.power(10, self.qs[self.lteff_name]) * u.K

    @property
    def lbol(self) -> Quantity["power"]:
        """Bolometric luminosities of the stellar track points."""
        return np.power(10, self.qs[self.llbol_name]) * u.Lsun

    @property
    def rad(self) -> Quantity["length"]:
        """Radii of the stellar track points."""
        return np.power(10, self.qs[self.lrad_name]) * u.Rsun

    @property
    def grav(self) -> Quantity["acceleration"]:
        """Surface gravities of the stellar track points."""
        return np.power(10, self.qs[self.lgrav_name]) * u.cm / u.s**2


@dataclass
class IsochroneSet:
    """
    Data class for storing a set of isochrones.

    Attributes:
        lages: Array of log_10(ages/yr) of the isochrones
        hdr_list: List of column headers for the isochrone data
        isos: List of Isochrone objects
    """

    lages: np.ndarray[np.float64]
    hdr_list: list[str]
    isos: list[Isochrone]
    # minimum and maximum masses that these isochrones can be applied to
    min_mass: Quantity["mass"]
    max_mass: Quantity["mass"]
    metallicity: float

    # names also stored as attributes for ease of access
    eep_name: str
    m0_name: str
    lteff_name: str
    llbol_name: str
    lrad_name: str
    lgrav_name: str
    elems: list[str]


@dataclass
class TrackSet:
    """
    Data class for storing a set of stellar evolution tracks.

    Attributes:
        masses: List of initial masses of the tracks
        hdr_list: List of column headers for the track data
        tracks: List of Track objects
        max_eep: Maximum equivalent evolutionary phase number across all tracks,
                 None for models without EEPs
    """

    masses: np.ndarray[Quantity["mass"]]
    min_ages: np.ndarray[Quantity["time"]]
    max_ages: np.ndarray[Quantity["time"]]
    hdr_list: list[str]
    tracks: list[StellarTrack]
    max_eep: int | None

    # names also stored as attributes for ease of access
    eep_name: str | None
    age_name: str
    lteff_name: str
    llbol_name: str
    lrad_name: str
    lgrav_name: str
    elems: list[str]
