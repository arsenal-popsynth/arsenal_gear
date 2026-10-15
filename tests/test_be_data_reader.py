"""
test_be_data_reader
===================

This file contains tests for the binary evolution data readers, run on small
fake BPASS and MPA/Bonn model trees (no downloads).
"""

# the tests read the converted tables through the readers' _load
# pylint: disable=protected-access

import io
import tarfile
import zipfile

import astropy.constants as const
import astropy.units as u
import numpy as np
import pandas as pd
import pytest

from arsenal_gear.binary_evolution import data_reader
from arsenal_gear.binary_evolution.be_data_structures import BinaryTrackSet
from arsenal_gear.binary_evolution.data_reader import (
    BPASSReader,
    MPAReader,
    _model_name,
    _parse_model_name,
)
from arsenal_gear.stellar_evolution.se_data_structures import TrackSet

N_BIN = 25  # points in each fake binary model
T_BIN = np.linspace(0, 1e7, N_BIN)  # binary models end at 10 Myr


def _write_bpass(path, time, p_mass, s_mass):
    """Write a fake BPASS model file (50 columns, see the BPASS manual)."""
    data = np.zeros((len(time), 50))
    data[:, 1], data[:, 2:5] = time, 1.5
    data[:, 5], data[:, 37] = p_mass, s_mass
    data[:, 46:49] = 0.7
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(path, data)


@pytest.fixture(name="bpass_dir")
def fixture_bpass_dir(tmp_path):
    """
    Fake BPASS tree at z014: single stars, plus binaries with a merger, a
    rejuvenated companion, a plain companion, and M1 = 2.3.
    """
    root = tmp_path / "BPASS"
    for m in (2.3, 3, 3.5, 5, 10, 20):
        time = np.linspace(0, 3e7, 40) + 1.3e5
        time[0] = 0
        _write_bpass(root / "NEWSINMODS/z014" / f"sneplot-z014-{m}", time, m, 0)

    bdir = root / "NEWBINMODS/NEWBINMODS/z014"
    winds = 20 - 0.01 * np.arange(N_BIN)
    gains = np.zeros(N_BIN)
    gains[[11, 12, 16]] = [1, 1, 2]  # primary gains mass at steps 10, 11, 15
    _write_bpass(bdir / "sneplot-z014-20-0.5-2", T_BIN, winds + np.cumsum(gains), 5.0)
    s_rej = np.full(N_BIN, 3.0)
    s_rej[5:] = [3.1, 3.25, 3.4] + [3.5] * (N_BIN - 8)  # accretes ~17%
    _write_bpass(bdir / "sneplot-z014-20-0.15-1", T_BIN, winds, s_rej)
    _write_bpass(bdir / "sneplot-z014-20-0.25-1", T_BIN, winds, 5.0)
    _write_bpass(bdir / "sneplot-z014-2.3-0.4-1.2", T_BIN, 2.3, 0.92)
    return root


@pytest.fixture(name="bpass")
def fixture_bpass(bpass_dir):
    """BPASS reader; constructing it converts both tables."""
    return BPASSReader(met=0.014, rootdir=bpass_dir)


def _model(reader, kind, name):
    """Rows of one model in a converted table."""
    table = reader._load(kind)
    return table[table.model == name].reset_index(drop=True)


class TestBPASSConversion:
    """Conversion of BPASS models into tables."""

    def test_tables(self, bpass):
        """Both tables have the shared schema and the expected model names."""
        singles = bpass._load("singles")
        binaries = bpass._load("binaries")
        star = ["mass", "logL", "logT", "logR"]
        assert list(singles.columns) == ["model", "time"] + star
        assert list(binaries.columns) == (
            ["model", "time"] + ["p_" + c for c in star] + ["s_" + c for c in star]
        )
        assert "00230" in set(singles.model)
        assert set(binaries.model) == {
            "00230_040_120",  # M1 = 2.3 rounds to 230, not 229
            "02000_015_100",
            "02000_025_100",
            "02000_050_200",
        }

    def test_merger(self, bpass):
        """The secondary is absorbed step by step, then no longer exists."""
        merger = _model(bpass, "binaries", "02000_050_200")
        net_gain = 1 - 0.01  # each merger step is net of the primary's wind loss
        expected = [5.0] * 11 + [5 - net_gain] + [5 - 2 * net_gain] * 4 + [0.0] * 9
        np.testing.assert_allclose(merger.s_mass[:N_BIN], expected)
        assert merger.loc[10:, ["s_logL", "s_logT", "s_logR"]].isna().all().all()
        assert len(merger) == N_BIN  # no single-star continuation after a merger

    @pytest.mark.parametrize("name", ["02000_015_100", "02000_025_100"])
    def test_continuation_moves_forward(self, bpass, name):
        """The companion's single-star continuation follows in system time."""
        system = _model(bpass, "binaries", name)
        assert len(system) > N_BIN
        assert np.all(np.diff(system.time) > 0)
        assert system.p_logL[N_BIN:].isna().all()



class TestConstruction:
    """Creating a reader produces the converted tables."""

    def test_builds_tables_and_keeps_raw_files(self, bpass_dir):
        """Both tables are written, and the raw model files stay in place."""
        reader = BPASSReader(met=0.014, rootdir=bpass_dir)
        for kind in ("singles", "binaries"):
            assert (reader.output_dir / f"{kind}_z014.pkl.gz").exists()
        assert len(list((bpass_dir / "NEWSINMODS/z014").iterdir())) == 6
        assert len(list((bpass_dir / "NEWBINMODS/NEWBINMODS/z014").iterdir())) == 4

    def test_existing_tables_reused(self, bpass_dir):
        """Existing tables are kept, unless overwrite=True rebuilds them."""
        reader = BPASSReader(met=0.014, rootdir=bpass_dir)
        paths = [reader.output_dir / f"{k}_z014.pkl.gz" for k in ("singles", "binaries")]
        mtimes = [path.stat().st_mtime_ns for path in paths]
        BPASSReader(met=0.014, rootdir=bpass_dir)
        assert [path.stat().st_mtime_ns for path in paths] == mtimes
        BPASSReader(met=0.014, rootdir=bpass_dir, overwrite=True)
        assert all(
            path.stat().st_mtime_ns != mtime for path, mtime in zip(paths, mtimes)
        )

    def test_convert_false(self, bpass_dir):
        """With convert=False, no tables are written."""
        reader = BPASSReader(met=0.014, rootdir=bpass_dir, convert=False)
        assert not list(reader.output_dir.iterdir())

    def test_no_model_files(self, tmp_path):
        """Model directories without model files are reported clearly."""
        for sub in ("NEWSINMODS/z014", "NEWBINMODS/NEWBINMODS/z014"):
            (tmp_path / sub).mkdir(parents=True)
        with pytest.raises(FileNotFoundError, match="NEWSINMODS"):
            BPASSReader(met=0.014, rootdir=tmp_path)


class TestBPASSGetData:
    """The download/unzip/untar cascade, with tiny archives."""

    @staticmethod
    def _release(root, models=True, tree=None):
        """
        A fake release zip holding the models tar file and an unneeded member.
        The tar holds the files of tree if given, otherwise a placeholder.
        """
        root.mkdir(parents=True, exist_ok=True)
        buf = io.BytesIO()
        with tarfile.open(fileobj=buf, mode="w:gz") as tar:
            if tree is not None:
                for sub in ("NEWSINMODS", "NEWBINMODS"):
                    tar.add(tree / sub, arcname=sub)
            else:
                member = "NEWBINMODS/NEWBINMODS/z014/x" if models else "OTHER/x"
                info = tarfile.TarInfo(member)
                info.size = 1
                tar.addfile(info, io.BytesIO(b"1"))
        with zipfile.ZipFile(root / BPASSReader.zip_name, "w") as zf:
            zf.writestr(BPASSReader.tar_name, buf.getvalue())
            zf.writestr("bpass_v2.2.1_imf135_300.tar.gz", b"unneeded")

    def test_zip_to_tables(self, tmp_path, bpass_dir):
        """
        From the release zip alone: only the models tar is extracted, the
        archives are deleted, the tables are built and the raw files are kept.
        """
        root = tmp_path / "release"
        self._release(root, tree=bpass_dir)
        reader = BPASSReader(met=0.014, rootdir=root)
        assert not (root / BPASSReader.zip_name).exists()
        assert not (root / BPASSReader.tar_name).exists()
        assert not (root / "bpass_v2.2.1_imf135_300.tar.gz").exists()
        assert len(list((root / "NEWBINMODS/NEWBINMODS/z014").iterdir())) == 4
        expected = BPASSReader(met=0.014, rootdir=bpass_dir)._load("binaries")
        pd.testing.assert_frame_equal(reader._load("binaries"), expected)

    def test_keep_archives(self, tmp_path):
        """delete_zip=False keeps the zip."""
        self._release(tmp_path)
        BPASSReader(rootdir=tmp_path, delete_zip=False, convert=False)
        assert (tmp_path / BPASSReader.zip_name).exists()
        assert (tmp_path / "NEWBINMODS").is_dir()

    def test_bad_layout(self, tmp_path):
        """A release without the models directory is reported clearly."""
        self._release(tmp_path, models=False)
        with pytest.raises(FileNotFoundError, match="layout"):
            BPASSReader(rootdir=tmp_path, convert=False)

    def test_no_download_by_default(self, tmp_path, monkeypatch):
        """Without local data, nothing is downloaded unless force_download=True."""

        def fail(*args, **kwargs):
            raise AssertionError("downloader should not be called")

        monkeypatch.setattr(data_reader, "downloader", fail)
        with pytest.raises(FileNotFoundError, match="force_download"):
            BPASSReader(rootdir=tmp_path)


def _write_mpa(path, n, m0):
    """Write a fake MPA/Bonn model with n points."""
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        {
            "star_age": np.linspace(0, 1e7, n),
            "star_mass": np.linspace(m0, 0.9 * m0, n),
            "log_L": np.linspace(3, 4, n),
            "log_Teff": np.linspace(4, 4.5, n),
            "log_R": np.linspace(0.5, 1, n),
        }
    ).to_pickle(path, compression="gzip")


@pytest.fixture(name="mpa")
def fixture_mpa(tmp_path):
    """MPA/Bonn reader for a fake MW tree, with the binaries converted."""
    root = tmp_path / "MPA"
    _write_mpa(root / "single_MW/1.000_0_compressed.pkl.gz", 30, 10.0)
    for name, n1, n2 in [
        ("1.300_0.500_1.000_compressed.pkl.gz", 10, 15),  # secondary outlives
        ("1.000_0.900_2.300_compressed.pkl.gz", 12, 8),  # merger
    ]:
        _write_mpa(root / "MW/primary/grid1" / name, n1, 20.0)
        _write_mpa(root / "MW/secondary/grid1" / name, n2, 10.0)
    return MPAReader(met="MW", rootdir=root)


class TestMPAConversion:
    """Conversion of MPA/Bonn models into tables."""

    def test_secondary_outlives_primary(self, mpa):
        """The primary keeps its final mass, and its other columns become NaN."""
        system = _model(mpa, "binaries", "01995_050_100")
        assert len(system) == 15
        np.testing.assert_allclose(system.p_mass[10:], 18.0)
        assert system.loc[10:, ["p_logL", "p_logT", "p_logR"]].isna().all().all()
        assert not system.filter(like="s_").isna().any().any()

    def test_merger(self, mpa):
        """After a merger, all of the secondary's columns are NaN."""
        system = _model(mpa, "binaries", "01000_090_230")  # logP = 2.3 rounds to 230
        assert len(system) == 12
        assert system.loc[8:, system.columns.str.startswith("s_")].isna().all().all()
        assert not system.filter(like="p_").isna().any().any()


def test_invalid_metallicity(tmp_path):
    """A metallicity the source doesn't provide is rejected."""
    with pytest.raises(ValueError):
        BPASSReader(met=0.5, rootdir=tmp_path)
    with pytest.raises(ValueError):
        MPAReader(met="M31", rootdir=tmp_path)


def test_load_before_convert(bpass_dir):
    """Loading a table that hasn't been converted names the method to run."""
    reader = BPASSReader(rootdir=bpass_dir, convert=False)
    with pytest.raises(FileNotFoundError, match="convert_single_data"):
        reader._load("singles")


@pytest.mark.parametrize(
    "params", [(20.0, None, None), (2.3, 0.4, 1.2), (100.0, 0.9, 3.8)]
)
def test_model_name_round_trip(params):
    """Model names decode back to the parameters they were built from."""
    assert _parse_model_name(_model_name(*params)) == params


class TestReadSingles:
    """Reading the singles table into a TrackSet."""

    def test_trackset(self, bpass):
        """One track per model, initial masses from the names, no EEPs."""
        tset = bpass.read_single_data()
        assert isinstance(tset, TrackSet)
        assert tset.max_eep is None and tset.eep_name is None
        np.testing.assert_allclose(
            tset.masses.to(u.Msun).value, [2.3, 3, 3.5, 5, 10, 20]
        )
        assert all(track.eeps is None for track in tset.tracks)

    def test_track_values(self, bpass):
        """Track properties match the table, and gravity follows from M and R."""
        track = bpass.read_single_data().tracks[-1]
        rows = _model(bpass, "singles", "02000")
        np.testing.assert_allclose(track.age.to(u.yr).value, rows.time)
        np.testing.assert_allclose(track.mass.to(u.Msun).value, rows.mass)
        np.testing.assert_allclose(track.lbol.to(u.Lsun).value, 10**rows.logL)
        g = const.G * track.mass / track.rad**2
        np.testing.assert_allclose(track.grav.to(u.cm / u.s**2).value,
                                   g.to(u.cm / u.s**2).value)

    def test_backward_time_steps_dropped(self, bpass_dir):
        """
        A row going back in time is dropped, while rows repeating a time are
        kept, so ages never decrease.
        """
        time = np.linspace(0, 3e7, 40)
        time[20] = time[19] * (1 - 1e-6)  # one tiny backward step
        time[30] = time[29]  # a repeated time
        _write_bpass(bpass_dir / "NEWSINMODS/z014/sneplot-z014-7", time, 7, 0)
        tset = BPASSReader(met=0.014, rootdir=bpass_dir).read_single_data()
        track = tset.tracks[list(tset.masses.value).index(7.0)]
        assert len(track.age) == 39
        assert np.all(np.diff(track.age.value) >= 0)


class TestReadBinaries:
    """Reading the binaries table into a BinaryTrackSet."""

    def test_initial_parameters(self, bpass):
        """Initial parameters are decoded from the model names."""
        bset = bpass.read_binary_data()
        assert isinstance(bset, BinaryTrackSet)
        binary = bset.tracks[list(bset.q0).index(0.5)]  # "02000_050_200"
        assert binary.primary.m0 == 20 * u.Msun
        assert binary.secondary.m0 == 10 * u.Msun
        assert binary.P0 == 100 * u.d
        assert binary.e0 == 0.0
        assert binary.a0.unit == u.au
        assert len(binary.primary.age) == len(binary.secondary.age)

    def test_merger(self, bpass):
        """After the merger, the secondary has no mass and no luminosity."""
        bset = bpass.read_binary_data()
        binary = bset.tracks[list(bset.q0).index(0.5)]
        assert binary.mass2[-1] == 0 * u.Msun
        assert np.isnan(binary.lbol2[10:].value).all()
        assert np.isnan(binary.grav2[-1].value)

    def test_mpa(self, mpa):
        """MPA/Bonn systems decode in the same way."""
        bset = mpa.read_binary_data()
        binary = bset.tracks[list(bset.q0).index(0.9)]  # "01000_090_230"
        assert binary.primary.m0 == 10 * u.Msun
        np.testing.assert_allclose(binary.P0.to(u.d).value, 10**2.3)
