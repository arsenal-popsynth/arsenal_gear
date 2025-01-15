#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 13 14:46:25 2025

@author: claude
"""

from typing import Type

import astropy.units as u
import numpy as np
from astropy.units import Quantity
from scipy.stats import rv_continuous, rv_histogram

import arsenal_gear

def test_fraction():
    # Create a StarPopulation
    mass = u.Msun*np.ones(100) 
    mass[50:] *= 20 # Change half the stars into massive stars
    metals = 0.1*np.ones(100)
    tform = u.Myr*np.zeros(100)
    stars = arsenal_gear.population.StarPopulation(mass=mass,metals=metals,tform=tform)
    #print(star)
    # Test the binary fraction sampler
    _binaries = arsenal_gear.dist_funcs.binaries.StepFraction([0, 1], 10*u.Msun, stars).sample()
    #print(_binaries)
    # Create the BinaryPopulation
    primaries = arsenal_gear.population.StarPopulation(mass=mass[_binaries],metals=metals[_binaries],tform=tform[_binaries])
    print(primaries)
    companions = arsenal_gear.population.StarPopulation(mass=0.5*mass[_binaries],metals=metals[_binaries],tform=tform[_binaries])
    semimajor = u.au*np.ones(len(primaries))
    eccentricity = np.zeros(len(primaries))
    binaries = arsenal_gear.population.BinaryPopulation(primary=primaries, secondary=companions, 
                                                        semimajor=semimajor, eccentricity=eccentricity)
    print(binaries)

test_fraction()