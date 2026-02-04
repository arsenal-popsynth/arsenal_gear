"""
test_isochrone
==========

This file contains unit tests for the isochrone interpolation methods
which are mostly contained in stellar_evolution/isochrone.py
"""

import astropy.units as u
import numpy as np
from numpy.testing import assert_array_less
from scipy.integrate import trapezoid

import arsenal_gear


def get_contiguous_regions(mask):
    """
    Return the start and end indices of contiguous True regions in a boolean mask.

    Helper function for identifying regions in initial stellar mass space that are
    missed by the isochrone interpolation, and are then separately integrated over
    to determine the fraction of luminosity that is missed by the interpolation.
    """
    dmask = np.diff(mask.astype(int))
    starts = np.where(dmask == 1)[0] + 1
    ends = np.where(dmask == -1)[0] + 1

    if mask[0]:
        starts = np.insert(starts, 0, 0)
    if mask[-1]:
        ends = np.append(ends, len(mask))

    return list(zip(starts, ends))


def integrate_mask(y, x, mask):
    """Integrate y over x where mask is True."""
    regions = get_contiguous_regions(mask)
    total_integral = 0.0
    for s, e in regions:
        if e - s > 1:
            total_integral += trapezoid(y[s:e], x[s:e])
    return total_integral


def test_mist_interp():
    """initialize stellar population with EEP interpolation
    and otherwise default parameters (which will specify MIST isochrones)
    test=True in isochrone interpolation to leave out nearest isochrone
    """
    sp = {
        "eep": arsenal_gear.StellarPopulation(interp_op="eep"),
        "iso": arsenal_gear.StellarPopulation(interp_op="iso", test=True),
    }

    lmissed = {"eep": [], "iso": []}
    L_err = {"eep": [], "iso": []}
    T_err = {"eep": [], "iso": []}
    T, L, nm, lum, lw_lerr = {}, {}, {}, {}, {}

    ais = np.arange(len(sp["iso"].iso.iset.lages))
    for ai in ais[2::6]:
        ai += 1
        t = (
            (1 + 1e-6)
            * np.power(10, np.array([sp["iso"].iso.iset.lages[ai]]) - 6)
            * u.Myr
        )

        ms = sp["iso"].iso.iset.isos[ai].qs["initial_mass"] * u.Msun
        xi = sp["iso"].imf.pdf(ms)

        L_ref = np.power(10, sp["iso"].iso.iset.isos[ai].qs["log_L"])
        T_ref = np.power(10, sp["iso"].iso.iset.isos[ai].qs["log_Teff"])
        lum_ref = trapezoid(L_ref * xi, ms.value)
        lw_teff_ref = trapezoid(L_ref * xi * T_ref, ms.value) / lum_ref

        for k in ["eep", "iso"]:
            T[k] = (sp[k].iso.teff(ms, t) / u.K).value
            L[k] = (sp[k].iso.lbol(ms, t) / u.Lsun).value

            # fraction of the total luminosity that
            # is missed by interpolation edge effects
            lmissed[k].append(integrate_mask(L_ref * xi, ms.value, L[k].mask) / lum_ref)

            nm[k] = np.logical_not(L[k].mask)

            lum[k] = integrate_mask(L[k] * xi, ms.value, nm[k])
            L_err[k].append(abs(lum[k] - lum_ref) / lum_ref)

            T_err[k].append(
                abs(
                    integrate_mask(L[k] * xi * T[k], ms.value, nm[k]) / lum[k]
                    - lw_teff_ref
                )
                / lw_teff_ref
            )

            lw_lerr[k] = (
                integrate_mask(np.abs(L[k] - L_ref) * xi, ms.value, nm[k]) / lum[k]
            )

    for k in ["eep", "iso"]:
        lmissed[k] = np.array(lmissed[k])
        L_err[k] = np.array(L_err[k])
        T_err[k] = np.array(T_err[k])
        for arr in lmissed[k], L_err[k], T_err[k]:
            # these don't seem like very stringent constraints, but there are regions
            # where things don't seem to work perfectly well in summary metrics, even
            # though the interpolation looks fine by eye
            assert_array_less(arr, 0.2)
            # the average errors over all time, especially if you weight by luminosity
            # over time, which we don't do, are quite small 3% or less
            assert np.average(arr) < 0.03


def test_lbol_methods_mist():
    """
    Test that compares both track-based and isochrone-based methods for
    interpolating the MIST isochronges, for both discrete and continuous
    forms of the stellar population.
    """
    int_ops = ["iso", "track"]
    discrete_ops = [True, False]
    outs = {}
    # array of times over which to compare bolometric luminosities
    tlin = np.logspace(5.1,9,20)*u.yr
    for int_op in int_ops:
        for discrete in discrete_ops:
            key = f"{int_op}_{'discrete' if discrete else 'continuous'}"
            sp = arsenal_gear.StellarPopulation(
                interp_op=int_op, discrete=discrete
            )
            outs[key] = sp.lbol(tlin)
    # compare all combinations of the four methods
    keys = list(outs.keys())
    n = len(keys)
    for i in range(n):
        for j in range(i + 1, n):
            rel_err = np.abs(1 - outs[keys[i]] / outs[keys[j]])
            # average relative error should be less than 4%
            assert np.average(rel_err) < 0.04