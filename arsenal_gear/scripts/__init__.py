"""
scripts
=========

Tools for running arsenal from the command line.
"""

import argparse

import astropy.units as u
import numpy as np

from arsenal_gear.formation.dist_funcs import *
from arsenal_gear.formation import *


def zams_single(
    IMF,
    *,
    N=None,
    mass=None,
    metals=0,
    min_mass=0.1 * u.Msun,
    max_mass=100 * u.Msun,
):
    """
    Create a zero-age main sequence single star population.
    """
    imf_ = IMF(min_mass, max_mass)
    if N is None and mass is None:
        raise ValueError("Either a total mass or number of stars must be specified.")
    if N is not None and mass is not None:
        raise ValueError("Fixing both total mass and number of stars is not possible.")
    if N is not None:  # We are operating in fixed-N mode
        masses = imf_.sample(N)
    else:  # We are operating in fixed-mass mode
        masses = imf_.sample_mass(mass)
    Mtot = np.sum(masses)
    return SinglePop(Mtot = Mtot,
                     Nstar = len(masses),
                     metallicity=metals * u.dimensionless_unscaled,
                     imf=imf_,
                     mmin=min_mass,
                     mmax=max_mass,
                     discrete=True)


def evolve_population(population, evolution, outputs, interval):
    """
    Evolve a population using a specified evolution model, recording specified
    outputs at specified intervals.
    """
    print(population, evolution, outputs, interval)


def main():
    """
    This is the primary method for running arsenal from the command line.
    """
    parser = argparse.ArgumentParser(
        description="Arsenal Gear is a Population Synthesis Tool"
    )
    parser.add_argument("paramfile")
    args = parser.parse_args()
    with open(args.paramfile, encoding="utf-8") as f:
        exec(f.read())
