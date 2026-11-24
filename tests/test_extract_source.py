"""
test_extract_source
===================

This file contains unit tests for the arsenal-extract-source command line
tool, which lives in arsenal_gear/scripts/extract_source.py
"""

# A test that takes a fixture names its parameter after that fixture, which
# pylint cannot tell apart from shadowing a global.
# pylint: disable=redefined-outer-name

import io
import tarfile

import h5py
import numpy as np
import pytest

import arsenal_gear
from arsenal_gear.file_io.metadata import construct_meta_dataset
from arsenal_gear.scripts.extract_source import extract_source, main


@pytest.fixture(scope="module")
def output_file(tmp_path_factory):
    """An output file with a full metadata group, shared across tests."""
    path = tmp_path_factory.mktemp("run") / "run.h5"
    with h5py.File(path, "w") as f:
        construct_meta_dataset(f, paramfile="mass = 1e4\n")
    return path


def make_tarball(members: dict) -> bytes:
    """Build a tar.gz in memory from a {name: bytes} dict."""
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        for name, data in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    return buf.getvalue()


def test_extracts_source_tree(output_file, tmp_path):
    dest = tmp_path / "src"
    assert main([str(output_file), "-o", str(dest)]) == 0

    extracted = dest / "arsenal_gear" / "file_io" / "metadata.py"
    original = f"{arsenal_gear.__path__[0]}/file_io/metadata.py"
    assert extracted.read_bytes() == open(original, "rb").read()


def test_default_output_dir(output_file, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert main([str(output_file)]) == 0
    assert (tmp_path / "run_source" / "arsenal_gear" / "__init__.py").is_file()


def test_info_does_not_extract(output_file, tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    assert main([str(output_file), "--info"]) == 0
    out = capsys.readouterr().out
    assert "rundate" in out
    assert "hostname" in out
    assert "git-commit" in out
    assert not list(tmp_path.iterdir())


def test_refuses_nonempty_output_dir(output_file, tmp_path, capsys):
    (tmp_path / "keep.txt").write_text("don't clobber me")
    assert main([str(output_file), "-o", str(tmp_path)]) == 1
    assert "--force" in capsys.readouterr().err
    assert not (tmp_path / "arsenal_gear").exists()


def test_force_extracts_into_nonempty_dir(output_file, tmp_path):
    (tmp_path / "keep.txt").write_text("don't clobber me")
    assert main([str(output_file), "-o", str(tmp_path), "--force"]) == 0
    assert (tmp_path / "arsenal_gear" / "__init__.py").is_file()
    assert (tmp_path / "keep.txt").read_text() == "don't clobber me"


def test_missing_source_dataset(tmp_path, capsys):
    path = tmp_path / "empty.h5"
    with h5py.File(path, "w") as f:
        f.create_dataset("metadata/rundate", data="2026-01-01T00:00:00")
    assert main([str(path), "-o", str(tmp_path / "src")]) == 1
    assert "metadata/source" in capsys.readouterr().err


def test_missing_file(tmp_path, capsys):
    assert main([str(tmp_path / "nope.h5")]) == 1
    assert "error" in capsys.readouterr().err


def test_unknown_archive_format(tmp_path, capsys):
    path = tmp_path / "zip.h5"
    with h5py.File(path, "w") as f:
        ds = f.create_dataset("metadata/source", data=np.void(b"PK\x03\x04"))
        ds.attrs["format"] = "zip"
    assert main([str(path), "-o", str(tmp_path / "src")]) == 1
    assert "zip" in capsys.readouterr().err


# Simulating a pre-filter Python on a modern one makes tarfile warn that it
# will filter by default in 3.14; real older Pythons don't emit this.
@pytest.mark.filterwarnings("ignore:Python 3.14 will:DeprecationWarning")
@pytest.mark.parametrize("has_filters", [True, False], ids=["filter", "fallback"])
def test_rejects_path_traversal(tmp_path, monkeypatch, has_filters):
    if not has_filters:  # exercise the pre-3.12 code path
        monkeypatch.delattr(tarfile, "data_filter", raising=False)
    blob = make_tarball({"../escaped.py": b"print('gotcha')\n"})
    with pytest.raises((tarfile.TarError, ValueError)):
        extract_source(blob, str(tmp_path / "src"))
    assert not (tmp_path / "escaped.py").exists()
