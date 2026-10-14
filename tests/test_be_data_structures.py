"""
test_be_data_structures
=======================

This file contains tests for the binary evolution data structures.
"""

import dataclasses

import astropy.units as u
import numpy as np
import pytest
from numpy.testing import assert_allclose

from arsenal_gear.binary_evolution.be_data_structures import BinaryTrack, BinaryTrackSet
from arsenal_gear.stellar_evolution.se_data_structures import StellarTrack, TrackSet

NAMES = {
    "age_name": "time",
    "lteff_name": "logT",
    "llbol_name": "logL",
    "lrad_name": "logR",
    "lgrav_name": "logg",
}


def _track(n, m0, mass=None, log_value=1.0):
    """A StellarTrack without EEPs, with n points and an optional orbit."""
    qs = {key: np.full(n, log_value) for key in ("logT", "logL", "logR", "logg")}
    qs["time"] = np.linspace(0, 1e7, n)
    qs["mass"] = np.full(n, m0) if mass is None else mass
    qs["period"] = np.linspace(10, 20, n)
    qs["ecc"] = np.linspace(0.1, 0.0, n)
    return StellarTrack(
        m0=m0 * u.Msun, eeps=None, mass_name="mass", elems=[], qs=qs, **NAMES
    )


def _binary(m1=20.0, q=0.5, period=10 * u.d, ecc=0.1, n=5, orbit_names=(None, None)):
    """A BinaryTrack whose stars share an n-point age grid."""
    return BinaryTrack(
        primary=_track(n, m1),
        secondary=_track(n, q * m1, log_value=2.0),
        q0=q,
        P0=period,
        e0=ecc,
        period_name=orbit_names[0],
        eccentricity_name=orbit_names[1],
    )


class TestBinaryTrack:
    """Tests for BinaryTrack."""

    def test_attributes(self):
        """The system's ages and initial parameters are stored as given."""
        binary = _binary()
        assert_allclose(binary.age.to(u.yr).value, binary.primary.age.to(u.yr).value)
        assert binary.secondary.m0 == binary.q0 * binary.primary.m0
        assert binary.P0 == 10 * u.d
        assert binary.e0 == 0.1

    def test_mismatched_age_grids(self):
        """Primary and secondary tracks of different lengths are rejected."""
        binary = _binary()
        with pytest.raises(ValueError):
            dataclasses.replace(binary, secondary=_track(4, 10.0))

    @pytest.mark.parametrize("prop", ["mass", "teff", "lbol", "rad", "grav"])
    def test_stellar_properties(self, prop):
        """Per-star properties come from the matching component track."""
        binary = _binary()
        for suffix, star in (("1", binary.primary), ("2", binary.secondary)):
            value = getattr(binary, prop + suffix)
            expected = getattr(star, prop)
            assert value.unit == expected.unit
            assert_allclose(value.value, expected.value)
        # the two stars hold different values, so they can't be mixed up
        assert not np.allclose(binary.lbol1.value, binary.lbol2.value)

    def test_mass_ratio(self):
        """The current mass ratio follows the current masses, e.g. during transfer."""
        binary = _binary()
        binary.secondary.qs["mass"] = np.linspace(10.0, 15.0, 5)
        assert_allclose(binary.mass_ratio, np.linspace(10.0, 15.0, 5) / 20.0)
        assert binary.q0 == 0.5

    def test_orbit_without_columns(self):
        """Current period, eccentricity and separation need their columns."""
        binary = _binary()
        for prop in ("period", "eccentricity", "semimajor_axis"):
            with pytest.raises(ValueError):
                getattr(binary, prop)

    def test_orbit_with_columns(self):
        """Current period (in days) and eccentricity are read from primary.qs."""
        binary = _binary(orbit_names=("period", "ecc"))
        assert_allclose(binary.period.to(u.d).value, np.linspace(10, 20, 5))
        assert_allclose(binary.eccentricity, np.linspace(0.1, 0.0, 5))

    def test_semimajor_axis_kepler(self):
        """The Sun and Earth on a one-year orbit are 1 AU apart."""
        n = 3
        earth = 3.003e-6
        binary = BinaryTrack(
            primary=_track(n, 1.0),
            secondary=_track(n, earth),
            q0=earth,
            P0=1 * u.yr,
            e0=0.0,
            period_name="period",
            eccentricity_name=None,
        )
        binary.primary.qs["period"] = np.full(n, (1 * u.yr).to(u.d).value)
        for a in (binary.semimajor_axis, binary.a0):
            assert a.unit == u.au
            assert_allclose(a.value, 1.0, rtol=1e-4)


def test_binary_trackset():
    """A BinaryTrackSet is a TrackSet whose per-system arrays line up."""
    systems = [_binary(m1=20.0, q=0.5), _binary(m1=40.0, q=0.9, ecc=0.0)]
    bset = BinaryTrackSet(
        masses=np.array([20.0, 40.0]) * u.Msun,
        min_ages=np.zeros(2) * u.yr,
        max_ages=np.full(2, 1e7) * u.yr,
        hdr_list=list(systems[0].primary.qs),
        tracks=systems,
        max_eep=None,
        eep_name=None,
        q0=np.array([b.q0 for b in systems]),
        P0=u.Quantity([b.P0 for b in systems]),
        e0=np.array([b.e0 for b in systems]),
        elems=[],
        **NAMES,
    )
    assert isinstance(bset, TrackSet)
    assert bset.max_eep is None and bset.eep_name is None
    for i, binary in enumerate(bset.tracks):
        assert binary.primary.m0 == bset.masses[i]
        assert binary.q0 == bset.q0[i]
        assert binary.P0 == bset.P0[i]
        assert binary.e0 == bset.e0[i]
    # the set's a0 matches each system's own a0
    assert bset.a0.unit == u.au
    assert_allclose(bset.a0.value, [b.a0.to(u.au).value for b in bset.tracks])
