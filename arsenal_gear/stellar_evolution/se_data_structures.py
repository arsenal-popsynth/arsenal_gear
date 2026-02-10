"""
se_data_structures.py
==========

This file defines various data structures used in the interpolation of 
results from stellar evolution models
"""

from dataclasses import dataclass
from astropy.units import Quantity
import numpy as np

@dataclass
class Isochrone:
    """
    Data class for storing isochrone data.
    
    Attributes:
        age: Age of the isochrone
        eep_name: variable name for equivalent evolutionary phase
        mini_name: variable name for initial stellar masses
        lteff_name: variable name for log_10(effective temperature)
        llbol_name: variable name for log_10(bolometric luminosities)
        qs: Dictionary of isochrone quantities. These should minimally include
            EEP, initial_mass, log_Teff, and log_L but can include others.
    """
    age: Quantity["time"]
    eep_name: str
    mini_name: str
    lteff_name: str
    llbol_name: str
    lrad_name: str
    lgrav_name: str
    # list of surface abundances, if available
    elems: list[str]
    qs: dict

# TODO(@ltancas): think of the best way to make age units consistent
#                 between the isochrone and stellar track data structures
@dataclass
class StellarTrack:
    """
    Data class for storing stellar track data.
    
    Attributes:
        mass: Initial mass of the stellar track
        eep_name: variable name for equivalent evolutionary phase
        age_name: variable name for age (in years by default)
        lteff_name: variable name for log_10(effective temperature)
        llbol_name: variable name for log_10(bolometric luminosities)
        qs: Dictionary of stellar track quantities. These should minimally include
            EEP, log_Teff, and log_L but can include others.
    """
    mass: Quantity["mass"]
    eeps: np.ndarray[int]
    age_name: str
    lteff_name: str
    llbol_name: str
    lrad_name: str
    lgrav_name: str
    # list of surface abundances, if available
    elems: list[str]
    qs: dict

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
    max_mass: Quantity["mass"]
    metallicity: float

    # names also stored as attributes for ease of access
    eep_name: str
    mini_name: str
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
        max_eep: Maximum equivalent evolutionary phase number across all tracks
    """
    masses: np.ndarray[Quantity["mass"]]
    min_ages: np.ndarray[Quantity["time"]]
    max_ages: np.ndarray[Quantity["time"]]
    hdr_list: list[str]
    tracks: list[StellarTrack]
    max_eep: int

    # names also stored as attributes for ease of access
    eep_name: str
    age_name: str
    lteff_name: str
    llbol_name: str
    lrad_name: str
    lgrav_name: str
    elems: list[str]
