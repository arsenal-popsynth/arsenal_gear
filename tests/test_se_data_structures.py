"""
test_se_data_structures
=======================

This file contains tests for the stellar evolution data structures:
their properties and support for models without EEPs on hand-built
objects, and the physical consistency and bookkeeping of the structures
built from (cached) MIST data.
"""

import astropy.constants as const
import astropy.units as u
import numpy as np
import pytest
from numpy.testing import assert_allclose

from arsenal_gear.stellar_evolution.data_reader import MISTReader
from arsenal_gear.stellar_evolution.se_data_structures import (
    Isochrone,
    StellarTrack,
    TrackSet,
)

NAMES = {
    "lteff_name": "logT",
    "llbol_name": "logL",
    "lrad_name": "logR",
    "lgrav_name": "logg",
}


def _qs():
    """Quantities for a three-point track or isochrone."""
    return {
        "age": np.array([0.0, 1e6, 2e6]),
        "mass": np.array([10.0, 9.5, 9.0]),
        "m0": np.array([1.0, 5.0, 10.0]),
        "EEP": np.array([1.0, 2.0, 3.0]),
        "logT": np.array([4.0, 4.1, 4.2]),
        "logL": np.array([3.0, 3.5, 4.0]),
        "logR": np.array([0.5, 0.6, 0.7]),
        "logg": np.array([4.0, 3.9, 3.8]),
    }


def _track(mass_name="mass"):
    """A StellarTrack without EEPs."""
    return StellarTrack(
        m0=10.0 * u.Msun,
        eeps=None,
        age_name="age",
        mass_name=mass_name,
        elems=[],
        qs=_qs(),
        **NAMES,
    )


def _isochrone(eep_name):
    """An Isochrone with or without EEPs."""
    return Isochrone(
        age=1 * u.Myr, eep_name=eep_name, m0_name="m0", elems=[], qs=_qs(), **NAMES
    )


def _check_log_properties(obj, qs):
    """teff, lbol, rad and grav undo the log10 and attach the right units."""
    assert_allclose(obj.teff.to(u.K).value, 10 ** qs["logT"])
    assert_allclose(obj.lbol.to(u.Lsun).value, 10 ** qs["logL"])
    assert_allclose(obj.rad.to(u.Rsun).value, 10 ** qs["logR"])
    assert_allclose(obj.grav.to(u.cm / u.s**2).value, 10 ** qs["logg"])


class TestStellarTrack:
    """Tests for StellarTrack without EEPs."""

    def test_properties(self):
        """Properties give the right values and units."""
        track = _track()
        qs = _qs()
        assert track.m0 == 10.0 * u.Msun
        assert_allclose(track.age.to(u.yr).value, qs["age"])
        assert_allclose(track.mass.to(u.Msun).value, qs["mass"])
        _check_log_properties(track, qs)

    def test_no_mass_column(self):
        """Asking for the current mass without a mass column raises."""
        with pytest.raises(ValueError):
            _ = _track(mass_name=None).mass


class TestIsochrone:
    """Tests for Isochrone with and without EEPs."""

    def test_properties(self):
        """Properties give the right values and units."""
        iso = _isochrone("EEP")
        qs = _qs()
        assert iso.age == 1 * u.Myr
        assert_allclose(iso.m0.to(u.Msun).value, qs["m0"])
        _check_log_properties(iso, qs)

    def test_properties_read_live_qs(self):
        """Properties read qs when called, so added or replaced columns are seen."""
        iso = _isochrone("EEP")
        iso.qs["logL"] = np.zeros(3)
        assert_allclose(iso.lbol.to(u.Lsun).value, 1.0)

    def test_eep(self):
        """With EEPs, eep returns integers."""
        eep = _isochrone("EEP").eep
        assert eep.dtype.kind == "i"
        assert list(eep) == [1, 2, 3]

    def test_no_eep(self):
        """Without EEPs, asking for them raises."""
        with pytest.raises(ValueError):
            _ = _isochrone(None).eep


def test_trackset_without_eeps():
    """A TrackSet can be built from tracks without EEPs."""
    tset = TrackSet(
        masses=np.array([10.0]) * u.Msun,
        min_ages=np.array([0.0]) * u.yr,
        max_ages=np.array([2e6]) * u.yr,
        hdr_list=list(_qs()),
        tracks=[_track()],
        max_eep=None,
        eep_name=None,
        age_name="age",
        elems=[],
        **NAMES,
    )
    assert tset.max_eep is None and tset.eep_name is None
    assert tset.tracks[0].eeps is None


# --- structures built from MIST data (cached by the isochrone tests) ---


@pytest.fixture(scope="module")
def mist_isochrones():
    """IsochroneSet for solar-metallicity, non-rotating MIST isochrones."""
    return MISTReader(met=0.0, interp_op="iso").read_iso_data()


@pytest.fixture(scope="module")
def mist_tracks():
    """TrackSet for solar-metallicity, non-rotating MIST tracks."""
    return MISTReader(met=0.0, interp_op="track").read_track_data()


def _assert_physically_consistent(obj, mass, grav_tol):
    """
    Check that lbol, teff and rad satisfy Stefan-Boltzmann, and that grav,
    mass and rad satisfy g = G M / R^2, to within the given tolerances in dex.
    These catch columns mapped to the wrong name field or wrong units.
    """
    l_sb = 4 * np.pi * obj.rad**2 * const.sigma_sb * obj.teff**4
    assert_allclose(np.log10((l_sb / obj.lbol).decompose().value), 0, atol=1e-3)
    g = const.G * mass / obj.rad**2
    assert_allclose(np.log10((g / obj.grav).decompose().value), 0, atol=grav_tol)


class TestMISTPhysics:
    """Physical consistency of the properties on MIST data."""

    def test_tracks(self, mist_tracks):
        """Track properties are consistent, using the current-mass property."""
        for track in mist_tracks.tracks:
            _assert_physically_consistent(track, track.mass, grav_tol=1e-3)

    def test_isochrones(self, mist_isochrones):
        """
        Isochrone properties are consistent. MIST interpolates each isochrone
        column separately, so log_g drifts from G M / R^2 by up to ~0.04 dex.
        """
        for iso in mist_isochrones.isos:
            mass = iso.qs["star_mass"] * u.Msun
            _assert_physically_consistent(iso, mass, grav_tol=0.05)


class TestMISTBookkeeping:
    """Consistency of the TrackSet and IsochroneSet built by the MIST reader."""

    def test_trackset(self, mist_tracks):
        """Set-level arrays line up with the tracks they describe."""
        tset = mist_tracks
        assert len(tset.tracks) == len(tset.masses)
        for i, track in enumerate(tset.tracks):
            age = track.age.to(u.yr).value
            assert track.m0 == tset.masses[i]
            assert_allclose(tset.min_ages[i].to(u.yr).value, age.min())
            assert_allclose(tset.max_ages[i].to(u.yr).value, age.max())
            # every column has one value per track point
            assert {len(col) for col in track.qs.values()} == {len(age)}
            # the interpolator stores each track in a max_eep-wide matrix
            assert len(age) <= tset.max_eep

    def test_track_eeps(self, mist_tracks):
        """Primary EEPs are increasing row numbers starting at the first row."""
        for track in mist_tracks.tracks:
            assert track.eeps[0] == 1
            assert np.all(np.diff(track.eeps) > 0)
        assert mist_tracks.max_eep == max(t.eeps[-1] for t in mist_tracks.tracks)

    def test_track_current_mass(self, mist_tracks):
        """Single stars start at their initial mass and only lose mass."""
        for track in mist_tracks.tracks:
            m = track.mass.to(u.Msun).value
            m0 = track.m0.to(u.Msun).value
            assert_allclose(m[0], m0, rtol=1e-3)
            assert np.all(m <= m0 * (1 + 1e-6))

    def test_isochroneset(self, mist_isochrones):
        """Set-level arrays line up with the isochrones they describe."""
        iset = mist_isochrones
        assert len(iset.isos) == len(iset.lages)
        assert iset.min_mass < iset.max_mass
        for lage, iso in zip(iset.lages, iset.isos):
            assert_allclose(iso.age.to(u.yr).value, 10**lage)
            n = len(iso.eep)
            assert {len(col) for col in iso.qs.values()} == {n}
            assert set(iset.elems) <= set(iso.qs)
