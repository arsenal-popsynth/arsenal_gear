"""
arsenal-gear
============

A lightweight population synthesis code with an emphasis on the quantities
relevant for stellar feedback from massive stars.
"""
import argparse
import astropy.units as u
import numpy as np
from astropy.units import Quantity

from . import dist_funcs, feedbacks, population, stellar_evolution

__version__ = '0.0.1'
__all__ = ['population', 'dist_funcs', 'feedbacks', 'stellar_evolution']

def check_arguments(args):
    if args.continuous and args.stochastic:
        raise ValueError("Cannot use both continuous and stochastic sampling simultaneously.")
    if not args.continuous and not args.stochastic:
        raise ValueError("Must specify either continuous or stochastic sampling.")
    if args.binary_fraction > 0:
        raise NotImplementedError("Binary systems are not yet implemented.")
    if args.continuous:
        raise NotImplementedError("Continuous sampling is not yet implemented.")
    return args

def main():
    parser = argparse.ArgumentParser(description="Arsenal Gear: A lightweight population synthesis code.")
    parser.add_argument("-s", "--stochastic", action="store_true", default=False, help="Use stochastic sampling of the IMF")
    parser.add_argument("-c", "--continuous", action="store_true", default=False, help="Use continuous sampling of the IMF")
    parser.add_argument("-v", "--version", action="version", version=__version__)
    parser.add_argument("-m", "--mass", type=float, default=1e6, help="Total mass of the stellar population in solar masses")
    parser.add_argument("-Z", "--metallicity", type=float, default=1, help="Metallicity of the stellar population in solar metallicity")
    parser.add_argument("-b", "--binary-fraction", type=float, default=0, help="Fraction of stars in binary systems")
    args = check_arguments(parser.parse_args())

    if args.stochastic:
        imf = dist_funcs.imf.Salpeter(min_mass=0.1*u.Msun, max_mass=100*u.Msun)
        masses = imf.sample_mass(args.mass * u.Msun)
        print(masses.max())
