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
class SingleStarTable:
    """
    Data class for storing evolutionary data for single stars evolved
    with binary evolution codes. 
    
    Attributes:
        masses: list of initial masses of the stars
        logA_name: variable name for log_10(age/yr)
        mass_name: variable name for stellar mass/Msun
        logT_name: variable name for log_10(effective temperature/K)
        logL_name: variable name for log_10(bolometric luminosities/Lsun)
        logR_name: variable name for log_10(stelalr radius/Rsun)
        qs: Dictionary of stellar track quantities. These should minimally include
            logA, mass, logT, logL, and logR but can include others.
    """
    masses: np.ndarray[Quantity["mass"]]
    logA_name: str
    mass_name: str
    logT_name: str
    logL_name: str
    logR_name: str
    qs: dict

@dataclass
class BinaryStarTable:
    """
    Data class for storing evolutionary data for binary stars evolved
    with binary evolution codes. 
    
    Attributes:
        masses: list of initial masses of the primaries
        mass_ratios: list of initial mass ratios companion mass/primary mass
        periods: list of initial orbital periods
        logA_name: variable name for log_10(age/yr)
        mass1_name: variable name for stellar mass/Msun for the primary
        logT1_name: variable name for log_10(effective temperature/K) for the primary
        logL1_name: variable name for log_10(bolometric luminosities/Lsun) for the primary
        logR1_name: variable name for log_10(stelalr radius/Rsun) for the primary
        mass2_name: variable name for stellar mass/Msun for the companion
        logT2_name: variable name for log_10(effective temperature/K) for the companion
        logL2_name: variable name for log_10(bolometric luminosities/Lsun) for the companion
        logR2_name: variable name for log_10(stelalr radius/Rsun) for the companion
        qs: Dictionary of stellar track quantities. These should minimally include
            logA, mass1, logT1, logL1, logR1, mass2, logT2, logL2, and logR2 
            but can include others.
    """
    masses: np.ndarray[Quantity["mass"]]
    mass_ratios: np.ndarray[np.float64]
    periods: np.ndarray[Quantity["time"]]
    logA_name: str
    mass1_name: str
    logT1_name: str
    logL1_name: str
    logR1_name: str
    mass2_name: str
    logT2_name: str
    logL2_name: str
    logR2_name: str
    qs: dict
