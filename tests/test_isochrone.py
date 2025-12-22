import astropy.units as u
import numpy as np
import pytest
from numpy.testing import assert_array_equal
from scipy.integrate import trapezoid

import arsenal_gear
from tqdm import tqdm

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
    for (s,e) in regions:
        if (e-s>1):
            total_integral += trapezoid(y[s:e], x[s:e])
    return total_integral

def test_mist_interp():
    # initialize stellar population with EEP interpolation
    # and otherwise default parameters (which will specify MIST isochrones)
    # test=Ture in isochrone interpolation to elave out nearest isochrone
    sp_eep = arsenal_gear.StellarPopulation(interp_op="eep")
    sp_iso = arsenal_gear.StellarPopulation(interp_op="iso",test=True)


    (lmissed_eep, lmissed_iso) = ([],[])
    (rel_err_eep, rel_err_iso) = ([],[])
    (L_err_eep, L_err_iso) = ([],[])
    (T_err_eep, T_err_iso) = ([],[])

    ais = np.arange(len(sp_iso.iso.iset.ages))
    for ai in ais[2::5]:
        ai += 1
        t = (1+1e-6)*np.power(10,np.array([sp_iso.iso.ages[ai]])-6)*u.Myr

        ms = sp_iso.iso.iset.isos[ai]["initial_mass"] * u.Msun
        xi = sp_iso.imf.pdf(ms)

        T_eep = (sp_eep.iso.teff(ms, t)/u.K).value
        L_eep = (sp_eep.iso.lbol(ms, t)/u.Lsun).value

        T_iso = (sp_iso.iso.teff(ms, t)/u.K).value
        L_iso = (sp_iso.iso.lbol(ms, t)/u.Lsun).value

        L_ref = np.power(10, sp_iso.iso.isos[ai]["log_L"])
        T_ref = np.power(10, sp_iso.iso.isos[ai]["log_Teff"])
        lum = trapezoid(L_ref*xi, ms.value)
        lw_teff = trapezoid(L_ref*xi*T_ref, ms.value)/lum

        # fraction of the total luminosity that
        # is missed by interpolation edge effects
        Lmiss_eep = integrate_mask(L_ref*xi, ms.value, L_eep.mask)/lum
        Lmiss_iso = integrate_mask(L_ref*xi, ms.value, L_iso.mask)/lum
        lmissed_eep.append(Lmiss_eep)
        lmissed_iso.append(Lmiss_iso)

        neepm = np.logical_not(L_eep.mask)
        nisom = np.logical_not(L_iso.mask)

        lum_eep = integrate_mask(L_eep*xi, ms.value, neepm)
        lum_iso = integrate_mask(L_iso*xi, ms.value, nisom)
        L_err_eep.append(abs(lum_eep-lum)/lum)
        L_err_iso.append(abs(lum_iso-lum)/lum)

        lw_teff_eep = integrate_mask(L_eep*xi*T_eep, ms.value, neepm)/lum_eep
        lw_teff_iso = integrate_mask(L_iso*xi*T_iso, ms.value, nisom)/lum_iso
        T_err_eep.append(abs(lw_teff_eep-lw_teff)/lw_teff)
        T_err_iso.append(abs(lw_teff_iso-lw_teff)/lw_teff)

        lw_lerr_eep = integrate_mask(np.abs(L_eep-L_ref)*xi, ms.value, neepm)/lum_eep
        lw_lerr_iso = integrate_mask(np.abs(L_iso-L_ref)*xi, ms.value, nisom)/lum_iso
        rel_err_eep.append(lw_lerr_eep)
        rel_err_iso.append(lw_lerr_iso)

    (lmissed_eep, lmissed_iso) = (np.array(lmissed_eep), np.array(lmissed_iso))
    (rel_err_eep, rel_err_iso) = (np.array(rel_err_eep), np.array(rel_err_iso))
    (L_err_eep, L_err_iso) = (np.array(L_err_eep), np.array(L_err_iso))
    (T_err_eep, T_err_iso) = (np.array(T_err_eep), np.array(T_err_iso))

    # these don't seem like very stringent constraints, but there are regions
    # where things don't seem to work perfectly well in summary metrics, even
    # though the interpolation looks fine by eye
    assert(np.all(lmissed_eep < 0.35))
    assert(np.all(lmissed_iso < 0.35))
    assert(np.all(L_err_eep   < 0.50))
    assert(np.all(T_err_eep   < 0.30))
    assert(np.all(T_err_iso   < 0.30))


    # the averages errors over all time, especially if you weight by luminosity
    # over time, which we don't do, are quite small 3% or less
    assert(np.average(lmissed_eep < 0.03))
    assert(np.average(lmissed_iso < 0.03))
    assert(np.average(L_err_eep)  < 0.03)
    assert(np.average(L_err_eep)  < 0.04)
    assert(np.average(T_err_eep)  < 0.03)
    assert(np.average(T_err_iso)  < 0.03)
