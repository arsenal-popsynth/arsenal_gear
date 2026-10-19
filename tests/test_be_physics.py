"""
test_be_physics
===============

Physics-based tests of the binary evolution data, run on the small BPASS
subset in tests/data/bpass (real BPASS v2.2.1 tracks at Z = 0.014): BPASS
single stars compared with MIST, and the physical consistency of BPASS
binaries after conversion (#44).
"""

from pathlib import Path

import astropy.units as u
import numpy as np
import pytest

from arsenal_gear.binary_evolution import BPASSReader
from arsenal_gear.stellar_evolution.data_reader import MISTReader
from arsenal_gear.utils.array_utils import make_monotonic_increasing

BPASS_SUBSET = Path(__file__).parent / "data" / "bpass"


@pytest.fixture(name="bpass", scope="module")
def fixture_bpass(tmp_path_factory):
    """BPASS reader for the committed subset, writing its tables to a tmp dir."""
    return BPASSReader(
        met=0.014, rootdir=BPASS_SUBSET, output_dir=tmp_path_factory.mktemp("bpass")
    )


@pytest.fixture(name="bpass_singles", scope="module")
def fixture_bpass_singles(bpass):
    """TrackSet of the subset's single stars."""
    return bpass.read_single_data()


@pytest.fixture(name="mist_tracks", scope="module")
def fixture_mist_tracks():
    """MIST tracks at [Fe/H] = 0 (Z = 0.0142), close to BPASS's Z = 0.014."""
    return MISTReader(met=0.0, interp_op="track").read_track_data()


def _track(tset, mass):
    """The track of a set with initial mass `mass` (Msun)."""
    return next(t for t in tset.tracks if np.isclose(t.m0.to_value(u.Msun), mass))


def _mist_from_zams(track):
    """MIST time since the ZAMS (EEP 202), log L and log Teff, and the TAMS age (EEP 454)."""
    zams, tams = track.eeps[1] - 1, track.eeps[3] - 1
    age = track.age.to_value(u.yr)
    time = age[zams:] - age[zams]
    log_l, log_t = (track.qs[k][zams:] for k in (track.llbol_name, track.lteff_name))
    return time, log_l, log_t, time[tams - zams]


def _bpass_ms_end(track):
    """
    End of the main sequence of a BPASS track (which has no EEPs): the point of
    fastest radius growth, i.e. the Hertzsprung gap, ignoring the track's end.
    """
    time = track.age.to_value(u.yr)
    growth = np.gradient(track.qs["logR"], time)
    return time[np.argmax(np.where(time < 0.98 * time[-1], growth, -np.inf))]


@pytest.mark.parametrize("mass", [1, 5, 20, 60])
class TestBPASSvsMIST:
    """
    BPASS and MIST single stars of the same mass evolve alike. The tolerances
    are about twice the differences measured on these tracks: on the main
    sequence |dlog L| <= 0.009 and |dlog Teff| <= 0.043 (BPASS slightly cooler,
    most at 60 Msun), and BPASS main sequences last 1.01-1.23 times MIST's.
    """

    def test_main_sequence(self, mass, bpass_singles, mist_tracks):
        """log L and log Teff agree at 20, 50 and 80% of the main sequence."""
        m_time, m_log_l, m_log_t, t_ms = _mist_from_zams(_track(mist_tracks, mass))
        bpass = _track(bpass_singles, mass)
        b_time = bpass.age.to_value(u.yr)
        ages = np.array([0.2, 0.5, 0.8]) * t_ms
        for b_q, m_q, tol in (
            (bpass.qs["logL"], m_log_l, 0.02),
            (bpass.qs["logT"], m_log_t, 0.08),
        ):
            diff = np.interp(ages, b_time, b_q) - np.interp(ages, m_time, m_q)
            assert np.all(np.abs(diff) < tol), diff

    def test_main_sequence_lifetime(self, mass, bpass_singles, mist_tracks):
        """BPASS's main sequence lasts 1.0-1.5 times MIST's."""
        t_ms = _mist_from_zams(_track(mist_tracks, mass))[3]
        ratio = _bpass_ms_end(_track(bpass_singles, mass)) / t_ms
        assert 1.0 <= ratio <= 1.5, ratio


@pytest.mark.parametrize("mass", [20, 60])
def test_total_lifetime(mass, bpass_singles, mist_tracks):
    """
    Massive stars live as long in both (measured: 1.06 at 20 Msun and 1.02 at
    60 Msun). Lower-mass BPASS tracks continue along the white dwarf cooling
    sequence, so their total ages aren't compared.
    """
    m_time = _mist_from_zams(_track(mist_tracks, mass))[0]
    b_time = _track(bpass_singles, mass).age.to_value(u.yr)
    assert 0.9 <= b_time[-1] / m_time[-1] <= 1.15


# --- binaries (#44) ---
# Checked on the full z014 grid (12,663 systems): time strictly increases, the
# first-row total mass is within 3% of M1 (1 + q0), after the primary's death its
# mass is held and the companion's properties stay defined, and the converter
# never adds mass: mergers stay within 0.04% of the initial M1 + M2. In 402
# low-mass mergers (M1 ~ 2-4 Msun) BPASS's own primary ends up more massive than
# the initial M1 + M2; those are flagged in BinaryTrackSet.mass_excess. Models
# already merged at t = 0 have a placeholder companion and no continuation.
# When the primary becomes a white dwarf (its final drop below log R = -1),
# BPASS no longer evolves the companion, which continues along its single star
# from then on: in the 1,775 such systems the companion's mass jumps by a median
# 1.6% at the join (> 10% in 347, from the single-star mass grid and mass
# transfer, > 50% in 7), and only 3 companions have already ended there.

BINARY_FILES = {
    "merger": "sneplot-z014-10-0.1-0",
    "rising_merger": "sneplot-z014-5-0.9-0.2",  # secondary ticks up while absorbed
    "excess_merger": "sneplot-z014-1.9-0.3-0",  # BPASS's primary exceeds M1 + M2
    "premerged": "sneplot-z014-1.1-0.5-0",  # merged at t = 0
    "plain": "sneplot-z014-10-0.1-0.8",
    "rejuvenated": "sneplot-z014-10-0.2-0.2",
    "white_dwarf": "sneplot-z014-3-0.7-3",  # the primary becomes a white dwarf
    "wd_dip": "sneplot-z014-2-0.3-0.8",  # ... after dipping below log R = -1
}


@pytest.fixture(name="bpass_binary_set", scope="module")
def fixture_bpass_binary_set(bpass):
    """BinaryTrackSet of the subset's binaries."""
    return bpass.read_binary_data()


def _index(bset, kind):
    """Position in the set of the system read from BINARY_FILES[kind]."""
    m1, q0, log_p = (float(x) for x in BINARY_FILES[kind].split("-")[-3:])
    return next(
        i
        for i, b in enumerate(bset.tracks)
        if np.isclose(b.primary.m0.to_value(u.Msun), m1)
        and np.isclose(b.q0, q0)
        and np.isclose(b.P0.to_value(u.d), 10**log_p)
    )


@pytest.fixture(name="bpass_binaries", scope="module")
def fixture_bpass_binaries(bpass_binary_set):
    """The subset's binaries, keyed by their kind."""
    return {kind: bpass_binary_set.tracks[_index(bpass_binary_set, kind)]
            for kind in BINARY_FILES}


def _raw(kind, cols):
    """Columns of the subset's raw BPASS file, by 0-indexed column number."""
    raw = BPASS_SUBSET / "NEWBINMODS/NEWBINMODS/z014" / BINARY_FILES[kind]
    return np.genfromtxt(raw, usecols=cols).T


def _bpass_total(kind):
    """M1 + M2 as written in the subset's raw BPASS file."""
    p_mass, s_mass = _raw(kind, (5, 37))
    return p_mass + s_mass


def _wd_row(kind):
    """Raw row at which the primary becomes a white dwarf (final log R < -1)."""
    log_r = _raw(kind, 2)
    return np.flatnonzero(log_r >= BPASSReader.wd_log_radius)[-1] + 1


@pytest.mark.parametrize("kind", list(BINARY_FILES))
def test_binary_basics(kind, bpass_binaries):
    """Time strictly increases, and the initial total mass is M1 (1 + q0)."""
    binary = bpass_binaries[kind]
    assert np.all(np.diff(binary.age.to_value(u.yr)) > 0)
    total = (binary.mass1 + binary.mass2)[0] / (binary.primary.m0 * (1 + binary.q0))
    assert abs(total.decompose().value - 1) < 0.05


def test_merger(bpass_binaries):
    """
    From the merger on, the secondary no longer shines and is absorbed until no
    mass is left, and the primary never exceeds the pre-merger total mass.
    """
    binary = bpass_binaries["merger"]
    gone = np.isnan(binary.lbol2.value)
    first = np.flatnonzero(gone)[0]
    assert first > 0 and np.all(gone[first:])
    assert binary.mass2[-1] == 0 * u.Msun
    before = binary.mass1[first - 1] + binary.mass2[first - 1]
    assert np.all(binary.mass1[first:] <= before * 1.001)  # rounding


@pytest.mark.parametrize("kind", ["merger", "rising_merger"])
def test_merger_mass_conservation(kind, bpass_binaries):
    """
    Mergers conserve mass: the total never exceeds the initial M1 + M2 by more
    than rounding (0.04% at most on the full grid), and the secondary keeps its
    BPASS mass until the merger starts, then has none and no longer shines.
    """
    binary = bpass_binaries[kind]
    m_init = binary.primary.m0 * (1 + binary.q0)
    total = (binary.mass1 + binary.mass2) / m_init
    assert np.all(total.decompose().value <= 1.001)
    gone = np.isnan(binary.lbol2.value)
    assert gone.any() and np.all(binary.mass2[gone] == 0)
    assert np.all(binary.mass2[~gone] > 0)


@pytest.mark.parametrize("kind", list(BINARY_FILES))
def test_no_mass_added(kind, bpass_binaries):
    """
    The converted M1 + M2 never exceeds BPASS's own, row by row, during the
    binary phase (up to a white dwarf primary's formation).
    """
    raw = _bpass_total(kind)
    if kind in ("white_dwarf", "wd_dip"):
        raw = raw[: _wd_row(kind) + 1]
    binary = bpass_binaries[kind]
    converted = (binary.mass1 + binary.mass2).to_value(u.Msun)[: len(raw)]
    assert np.all(converted <= raw * (1 + 1e-9))


@pytest.mark.parametrize("kind", ["plain", "rejuvenated"])
def test_companion_continuation(kind, bpass_binaries):
    """
    After the primary's death its mass is held and it no longer shines, while
    the companion continues with defined properties, later in time.
    """
    binary = bpass_binaries[kind]
    dead = np.isnan(binary.lbol1.value)
    join = np.flatnonzero(dead)[0]
    assert join > 0 and np.all(dead[join:])
    assert np.all(binary.mass1[join:] == binary.mass1[join - 1])
    for prop in (binary.lbol2, binary.teff2, binary.rad2):
        assert np.all(np.isfinite(prop[join:].value))


def test_rejuvenated_companion(bpass_binaries, bpass_singles):
    """
    A rejuvenated companion continues along the single star of its maximum
    mass, as if born when it reached that mass: its ages and masses after the
    join are those of that single star, shifted by the birth time.
    """
    binary = bpass_binaries["rejuvenated"]
    time = binary.age.to_value(u.yr)
    mass2 = binary.mass2.to_value(u.Msun)
    join = np.flatnonzero(np.isnan(binary.lbol1.value))[0]
    m2_eff = mass2[:join].max()
    t_birth = time[np.flatnonzero(mass2[:join] == m2_eff)[0]]
    single = min(bpass_singles.tracks, key=lambda t: abs(t.m0.to_value(u.Msun) - m2_eff))
    s_time = single.age.to_value(u.yr)
    after = s_time > time[join - 1] - t_birth
    np.testing.assert_allclose(time[join:], t_birth + s_time[after])
    np.testing.assert_allclose(mass2[join:], single.mass.to_value(u.Msun)[after])


@pytest.mark.parametrize("kind", ["white_dwarf", "wd_dip"])
def test_white_dwarf_join(kind, bpass_binaries, bpass_singles):
    """
    When the primary becomes a white dwarf, the binary phase ends: the
    companion continues along its single star (as if born at t = 0, without
    rejuvenation) from the white dwarf's formation, with little jump in mass,
    and the primary keeps its whole raw track, cooling sequence included.
    """
    binary = bpass_binaries[kind]
    time, log_r1, log_l1, p_mass, m2_raw = _raw(kind, (1, 2, 4, 5, 37))
    i_wd = _wd_row(kind)
    age = binary.age.to_value(u.yr)
    mass2 = binary.mass2.to_value(u.Msun)
    # the binary phase is the raw rows up to the white dwarf's formation (with
    # tied times spread, by less than the raw times' 7 significant digits)
    np.testing.assert_allclose(age[: i_wd + 1], time[: i_wd + 1], rtol=1e-6)
    np.testing.assert_allclose(mass2[: i_wd + 1], m2_raw[: i_wd + 1])
    # every raw row of the primary is kept: the binary phase, then the cooling
    # rows at their (strictly increasing) times
    cooling = make_monotonic_increasing(
        np.arange(len(time) - i_wd, dtype=float), time[i_wd:].copy()
    )
    rows = np.isin(age, cooling[1:]) & (age > age[i_wd])
    rows[: i_wd + 1] = True
    assert rows.sum() == len(time)
    np.testing.assert_allclose(binary.mass1[rows].to_value(u.Msun), p_mass)
    np.testing.assert_allclose(np.log10(binary.lbol1[rows].to_value(u.Lsun)), log_l1)
    np.testing.assert_allclose(np.log10(binary.rad1[rows].to_value(u.Rsun)), log_r1)
    # then the companion follows its single star, shifted by nothing
    single = min(bpass_singles.tracks, key=lambda t: abs(t.m0.to_value(u.Msun) - m2_raw[0]))
    s_time = single.age.to_value(u.yr)
    after = age > time[i_wd]
    on_track = np.isin(age[after], s_time)
    assert on_track.sum() == np.sum(s_time > time[i_wd])
    np.testing.assert_allclose(mass2[after][on_track], single.mass.to_value(u.Msun)[s_time > time[i_wd]])
    # with a small jump in mass at the join (measured: 0% and 3.5%)
    assert abs(mass2[i_wd + 1] / mass2[i_wd] - 1) < 0.07


def test_white_dwarf_companion_evolves(bpass_binaries):
    """
    The companion of 3-0.7-3 now goes through its own giant phase and
    becomes a white dwarf (at ~1.3e9 yr), long before the binary model ends
    (~1e11 yr), while the first white dwarf is still cooling.
    """
    binary = bpass_binaries["white_dwarf"]
    age = binary.age.to_value(u.yr)
    log_r2 = np.log10(binary.rad2.to_value(u.Rsun))
    end = _raw("white_dwarf", 1)[-1]
    giant = age[log_r2 > 1]
    assert giant.size and giant[0] < end
    second_wd = age[(age > giant[0]) & (log_r2 < BPASSReader.wd_log_radius)][0]
    assert second_wd < 0.1 * end
    assert np.isfinite(binary.lbol1[age < end].value).all()


def test_mass_excess_flag(bpass_binary_set):
    """Only the merger whose BPASS primary outgrows the initial M1 + M2 is flagged."""
    flagged = {
        kind for kind in BINARY_FILES
        if bpass_binary_set.mass_excess[_index(bpass_binary_set, kind)]
    }
    assert flagged == {"excess_merger"}
    assert len(bpass_binary_set.mass_excess) == len(bpass_binary_set.tracks)


def test_already_merged(bpass_binaries):
    """
    A model merged at t = 0 has no companion at any time (BPASS's placeholder
    isn't treated as a star) and no single-star continuation.
    """
    binary = bpass_binaries["premerged"]
    assert np.all(binary.mass2 == 0 * u.Msun)
    assert np.all(np.isnan(binary.lbol2.value))
    raw_rows = len(np.genfromtxt(
        BPASS_SUBSET / "NEWBINMODS/NEWBINMODS/z014" / BINARY_FILES["premerged"], usecols=1
    ))
    assert len(binary.age) == raw_rows
    assert np.all(np.isfinite(binary.lbol1.value))
