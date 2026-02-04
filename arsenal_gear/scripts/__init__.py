# pylint: disable=exec-used,undefined-variable
"""
scripts
=========

Tools for running arsenal from the command line.
"""

import argparse

import astropy.units as u
import numpy as np

from arsenal_gear import *


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
        exec(f.read(), globals())
        for t0, t1 in zip(out_times[:-1], out_times[1:]):
            print(t0.to("Myr"), t1.to("Myr"), *[qty(t0, t1) for qty in out_qtys])
