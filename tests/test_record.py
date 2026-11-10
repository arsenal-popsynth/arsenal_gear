"""
test_record
==========

This file contains unit tests for writing timestepped records into hdf5
outputs, which lives in arsenal_gear/file_io/record.py
"""

# A test that takes a fixture names its parameter after that fixture, which
# pylint cannot tell apart from shadowing a global.
# pylint: disable=redefined-outer-name

import astropy.units as u
import h5py
import numpy as np
import pytest
from numpy.testing import assert_allclose

from arsenal_gear.file_io.record import add_record, write_recordset


@pytest.fixture
def h5file(tmp_path):
    """A fresh, writable hdf5 file."""
    with h5py.File(tmp_path / "test.h5", "w") as f:
        yield f


@pytest.fixture
def times():
    return np.logspace(-3, 2, 10) * u.Myr


@pytest.fixture
def energies():
    return np.linspace(1, 10, 10) * 1e51 * u.erg


# --- time handling --------------------------------------------------------


@pytest.mark.parametrize("bad_time", [np.arange(5.0), [1.0, 2.0], 1.0, None], ids=repr)
def test_rejects_other_time_types(h5file, energies, bad_time):
    with pytest.raises(TypeError):
        add_record(h5file, "a", bad_time, energies)
    assert "data" not in h5file


def test_rejects_incompatible_time_units(h5file, times, energies):
    with pytest.raises(u.UnitConversionError):
        add_record(h5file, "a", times, energies, time_units=u.cm)


def test_link_makes_softlink(h5file, times, energies):
    add_record(h5file, "a", times, energies)
    add_record(h5file, "b", "a", energies)

    link = h5file["data/b"].get("time", getlink=True)
    assert isinstance(link, h5py.SoftLink)
    assert link.path == "/data/a/time"
    assert_allclose(h5file["data/b/time"][()], times.value)
    assert h5file["data/b/time"].attrs["units"] == h5file["data/a/time"].attrs["units"]


def test_link_written_before_target(h5file, times, energies):
    """A soft link may dangle until the record it points to is written."""
    add_record(h5file, "b", "a", energies)
    with pytest.raises(KeyError):
        h5file["data/b/time"]  # pylint: disable=pointless-statement
    add_record(h5file, "a", times, energies)
    assert_allclose(h5file["data/b/time"][()], times.value)


def test_link_chain(h5file, times, energies):
    add_record(h5file, "a", times, energies)
    add_record(h5file, "b", "a", energies)
    add_record(h5file, "c", "b", energies)
    assert_allclose(h5file["data/c/time"][()], times.value)


def test_link_ignores_time_units(h5file, times, energies):
    """A linked record shares its target's times, units and all."""
    add_record(h5file, "a", times, energies)
    add_record(h5file, "b", "a", energies, time_units=u.yr)
    assert u.Unit(h5file["data/b/time"].attrs["units"]) == u.Myr


def test_same_name_twice_fails(h5file, times, energies):
    add_record(h5file, "a", times, energies)
    with pytest.raises(ValueError):
        add_record(h5file, "a", times, energies)


# --- add_record -----------------------------------------------------------


def test_add_record_keeps_units(h5file, times, energies):
    add_record(h5file, "E", times, energies)
    qty = h5file["data/E/qty"]
    assert_allclose(qty[()], energies.value)
    assert u.Unit(qty.attrs["units"]) == u.erg
    assert_allclose(h5file["data/E/time"][()], times.value)
    assert u.Unit(h5file["data/E/time"].attrs["units"]) == u.Myr


def test_add_record_converts_units(h5file, times, energies):
    add_record(
        h5file,
        "E",
        times,
        energies,
        record_units=u.J,
        time_units=u.Gyr,
    )
    qty = h5file["data/E/qty"]
    assert_allclose(qty[()], energies.to_value(u.J))
    assert u.Unit(qty.attrs["units"]) == u.J
    assert_allclose(h5file["data/E/time"][()], times.to_value(u.Gyr))
    assert u.Unit(h5file["data/E/time"].attrs["units"]) == u.Gyr


def test_add_record_roundtrips_to_quantity(h5file, times, energies):
    """What's written can be turned back into the original quantity."""
    add_record(h5file, "E", times, energies, record_units=u.J)
    qty = h5file["data/E/qty"]
    restored = qty[()] * u.Unit(qty.attrs["units"])
    assert_allclose(restored.to_value(u.erg), energies.to_value(u.erg))


def test_add_record_with_linked_time(h5file, times, energies):
    masses = np.ones(10) * u.Msun
    add_record(h5file, "E", times, energies)
    add_record(h5file, "m", "E", masses)

    assert_allclose(h5file["data/m/qty"][()], masses.value)
    assert u.Unit(h5file["data/m/qty"].attrs["units"]) == u.Msun
    link = h5file["data/m"].get("time", getlink=True)
    assert isinstance(link, h5py.SoftLink)
    assert link.path == "/data/E/time"


def test_add_record_survives_reopen(tmp_path, times, energies):
    path = tmp_path / "test.h5"
    with h5py.File(path, "w") as f:
        add_record(f, "E", times, energies)
        add_record(f, "E2", "E", energies)
    with h5py.File(path, "r") as f:
        assert_allclose(f["data/E/qty"][()], energies.value)
        assert_allclose(f["data/E2/time"][()], times.value)


# --- write_recordset ------------------------------------------------------


def is_link(f, name):
    return isinstance(f[f"data/{name}"].get("time", getlink=True), h5py.SoftLink)


def test_recordset_links_duplicate_times(h5file, times, energies):
    masses = np.ones(10) * u.Msun
    other_times = np.linspace(0, 1, 10) * u.Myr
    records = {
        "E": {"data": energies, "time": times},
        "m": {"data": masses, "time": times.copy()},
        "E2": {"data": energies, "time": other_times},
        "m2": {"data": masses, "time": other_times.copy()},
    }
    write_recordset(h5file, records)

    assert not is_link(h5file, "E")
    assert not is_link(h5file, "E2")
    assert h5file["data/m"].get("time", getlink=True).path == "/data/E/time"
    assert h5file["data/m2"].get("time", getlink=True).path == "/data/E2/time"
    assert_allclose(h5file["data/m/time"][()], times.value)
    assert_allclose(h5file["data/m2/time"][()], other_times.value)
    assert_allclose(h5file["data/m/qty"][()], masses.value)
    assert u.Unit(h5file["data/m/qty"].attrs["units"]) == u.Msun


def test_recordset_does_not_modify_input(h5file, times, energies):
    records = {
        "a": {"data": energies, "time": times},
        "b": {"data": energies, "time": times},
    }
    write_recordset(h5file, records)
    assert records["b"]["time"] is times


def test_recordset_passes_units(h5file, times, energies):
    records = {
        "E": {
            "data": energies,
            "time": times,
            "record_units": u.J,
            "time_units": u.yr,
        }
    }
    write_recordset(h5file, records)
    assert_allclose(h5file["data/E/qty"][()], energies.to_value(u.J))
    assert u.Unit(h5file["data/E/qty"].attrs["units"]) == u.J
    assert_allclose(h5file["data/E/time"][()], times.to_value(u.yr))
    assert u.Unit(h5file["data/E/time"].attrs["units"]) == u.yr


def test_recordset_keeps_differing_time_units_separate(h5file, times, energies):
    """Linking would silently drop the second record's requested time units."""
    records = {
        "a": {"data": energies, "time": times},
        "b": {"data": energies, "time": times, "time_units": u.yr},
    }
    write_recordset(h5file, records)
    assert not is_link(h5file, "b")
    assert u.Unit(h5file["data/b/time"].attrs["units"]) == u.yr


def test_recordset_accepts_explicit_links(h5file, times, energies):
    records = {
        "b": {"data": energies, "time": "a"},
        "a": {"data": energies, "time": times},
    }
    write_recordset(h5file, records)
    assert is_link(h5file, "b")
    assert_allclose(h5file["data/b/time"][()], times.value)


def test_recordset_rejects_bad_time(h5file, energies):
    with pytest.raises(TypeError):
        write_recordset(h5file, {"a": {"data": energies, "time": [1.0, 2.0]}})
