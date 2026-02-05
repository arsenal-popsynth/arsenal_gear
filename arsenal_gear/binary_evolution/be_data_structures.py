"""
be_data_structures.py
==========

This file defines the data structures used to access the stellar
properties derived from binary evolution models.
"""

from dataclasses import dataclass
from astropy.units import Quantity
import numpy as np


@dataclass
class SingleStarTrack:
    """
    Data class for storing evolutionary data for single stars evolved
    with binary evolution codes. Note that this is not the same as the
    StellarTrack class implemented in se_data_structures.py
    
    Attributes:
        model_name: model name, based on the mass of the star
        age_name: variable name for age (in years by default)
        lteff_name: variable name for log_10(effective temperature/K)
        llbol_name: variable name for log_10(bolometric luminosities/Lsun)
        qs: Dictionary of stellar track quantities. These should minimally include
            log_Teff and log_L but can include others.
    """
    mass: Quantity["mass"] # convert from mass?
    age_name: str
    lteff_name: str
    llbol_name: str
    qs: dict

@dataclass
class BinaryStarTrack:
    """
    Data class for storing evolutionary data for binary stars evolved
    with binary evolution codes. 
    
    Attributes:
        model_name: model name, based on primary mass, period, and mass ratio
        lteff_name: variable name for log_10(effective temperature/K)
        llbol_name: variable name for log_10(bolometric luminosities/Lsun)
        qs: Dictionary of stellar track quantities. These should minimally include
            log_Teff and log_L but can include others.
    """
    mass: Quantity["mass"]
    age_name: str
    lteff_name: str
    llbol_name: str
    qs: dict


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
