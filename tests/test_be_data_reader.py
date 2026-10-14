"""
test_be_data_reader
===================

This file contains tests for the binary evolution data readers, run on small
fake BPASS and MPA/Bonn model trees (no downloads).
"""

import io
import tarfile
import zipfile

import numpy as np
import pandas as pd
import pytest

from arsenal_gear.binary_evolution import data_reader
from arsenal_gear.binary_evolution.data_reader import BPASSReader, MPAReader

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
    """BPASS reader with both tables converted."""
    reader = BPASSReader(met=0.014, rootdir=bpass_dir)
    reader.convert_single_data()
    reader.convert_binary_data()
    return reader


def _model(reader, kind, name):
    """Rows of one model in a converted table."""
    table = reader._load(kind)  # pylint: disable=protected-access
    return table[table.model == name].reset_index(drop=True)


class TestBPASSConversion:
    """Conversion of BPASS models into tables."""

    def test_tables(self, bpass):
        """Both tables have the shared schema and the expected model names."""
        singles = bpass._load("singles")  # pylint: disable=protected-access
        binaries = bpass._load("binaries")  # pylint: disable=protected-access
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

    def test_overwrite(self, bpass_dir):
        """Existing tables are only replaced with overwrite=True."""
        reader = BPASSReader(met=0.014, rootdir=bpass_dir)
        reader.convert_single_data()
        path = reader.output_dir / "singles_z014.pkl.gz"
        mtime = path.stat().st_mtime_ns
        reader.convert_single_data()
        assert path.stat().st_mtime_ns == mtime
        BPASSReader(met=0.014, rootdir=bpass_dir, overwrite=True).convert_single_data()
        assert path.stat().st_mtime_ns != mtime


class TestBPASSGetData:
    """The download/unzip/untar cascade, with tiny archives."""

    @staticmethod
    def _release(root, models=True):
        """A fake release zip holding the models tar file and an unneeded member."""
        root.mkdir(parents=True, exist_ok=True)
        buf = io.BytesIO()
        with tarfile.open(fileobj=buf, mode="w:gz") as tar:
            member = "NEWBINMODS/NEWBINMODS/z014/x" if models else "OTHER/x"
            info = tarfile.TarInfo(member)
            info.size = 1
            tar.addfile(info, io.BytesIO(b"1"))
        with zipfile.ZipFile(root / BPASSReader.zip_name, "w") as zf:
            zf.writestr(BPASSReader.tar_name, buf.getvalue())
            zf.writestr("bpass_v2.2.1_imf135_300.tar.gz", b"unneeded")

    def test_zip_to_models(self, tmp_path):
        """Only the models tar is extracted, and the archives are deleted."""
        self._release(tmp_path)
        BPASSReader(rootdir=tmp_path)
        assert (tmp_path / "NEWBINMODS").is_dir()
        assert not (tmp_path / BPASSReader.zip_name).exists()
        assert not (tmp_path / BPASSReader.tar_name).exists()
        assert not (tmp_path / "bpass_v2.2.1_imf135_300.tar.gz").exists()

    def test_keep_archives(self, tmp_path):
        """delete_zip=False keeps the zip."""
        self._release(tmp_path)
        BPASSReader(rootdir=tmp_path, delete_zip=False)
        assert (tmp_path / BPASSReader.zip_name).exists()

    def test_bad_layout(self, tmp_path):
        """A release without the models directory is reported clearly."""
        self._release(tmp_path, models=False)
        with pytest.raises(FileNotFoundError):
            BPASSReader(rootdir=tmp_path)

    def test_no_download_by_default(self, tmp_path, monkeypatch):
        """Without local data, nothing is downloaded unless force_download=True."""

        def fail(*args, **kwargs):
            raise AssertionError("downloader should not be called")

        monkeypatch.setattr(data_reader, "downloader", fail)
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
    reader = MPAReader(met="MW", rootdir=root)
    reader.convert_binary_data()
    return reader


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


def test_load_before_convert(tmp_path):
    """Loading a table that hasn't been converted names the method to run."""
    reader = BPASSReader(rootdir=tmp_path)
    with pytest.raises(FileNotFoundError, match="convert_single_data"):
        reader._load("singles")  # pylint: disable=protected-access
