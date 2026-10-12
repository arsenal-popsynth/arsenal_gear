"""
test_file_io
============

This file contains tests for the archive handling utilities in utils.file_io.
"""

import tarfile
from zipfile import ZipFile

import pytest

from arsenal_gear.utils import is_valid_tar, untar, unzip

MEMBERS = {"a.txt": "alpha", "b.txt": "beta"}


def _write_members(directory):
    """Write the test member files to directory and return their paths."""
    directory.mkdir(parents=True, exist_ok=True)
    paths = []
    for name, content in MEMBERS.items():
        path = directory / name
        path.write_text(content)
        paths.append(path)
    return paths


TAR_MODES = {"tar.gz": "w:gz", "txz": "w:xz", "tar": "w"}


def _make_tar(directory, mode):
    """Build a tar archive containing MEMBERS in directory/archive."""
    paths = _write_members(directory / "src")
    archive = directory / "archive" / "data.tar"
    archive.parent.mkdir()
    with tarfile.open(archive, mode) as tar:
        for path in paths:
            tar.add(path, arcname=path.name)
    return archive


def _make_zip(directory):
    """Build a zip archive containing MEMBERS in directory/archive."""
    paths = _write_members(directory / "src")
    archive = directory / "archive" / "data.zip"
    archive.parent.mkdir()
    with ZipFile(archive, "w") as zf:
        for path in paths:
            zf.write(path, arcname=path.name)
    return archive


@pytest.fixture(params=list(TAR_MODES) + ["zip"])
def archive(tmp_path, request):
    """(archive path, extraction function, delete kwarg) for each archive type."""
    if request.param == "zip":
        return _make_zip(tmp_path), unzip, "delete_zip"
    return _make_tar(tmp_path, TAR_MODES[request.param]), untar, "delete_tar"


class TestExtract:
    """Tests shared by untar and unzip."""

    @pytest.mark.parametrize("as_str", [False, True], ids=["Path", "str"])
    def test_extract_all(self, archive, as_str, tmp_path):
        """All members are extracted, whether paths are given as str or Path."""
        (path, extract, _) = archive
        out = tmp_path / "out"
        if as_str:
            res = extract(str(path), str(out))
        else:
            res = extract(path, out)
        assert res == out
        for name, content in MEMBERS.items():
            assert (out / name).read_text() == content
        assert path.exists()

    def test_default_target_dir(self, archive):
        """With no target_dir, members land next to the archive."""
        (path, extract, _) = archive
        res = extract(path)
        assert res == path.parent
        assert all((path.parent / name).exists() for name in MEMBERS)

    def test_target_file(self, archive, tmp_path):
        """Only the requested member is extracted."""
        (path, extract, _) = archive
        out = tmp_path / "out"
        extract(path, out, target_file="a.txt")
        assert (out / "a.txt").exists()
        assert not (out / "b.txt").exists()

    def test_missing_target_file(self, archive, tmp_path):
        """Requesting a member that isn't in the archive raises."""
        (path, extract, _) = archive
        with pytest.raises(FileNotFoundError):
            extract(path, tmp_path / "out", target_file="missing.txt")

    def test_inspect(self, archive, tmp_path, capsys):
        """inspect prints member names and extracts nothing."""
        (path, extract, _) = archive
        out = tmp_path / "out"
        extract(path, out, inspect=True)
        printed = capsys.readouterr().out
        assert all(name in printed for name in MEMBERS)
        assert not out.exists()

    def test_delete(self, archive, tmp_path):
        """The archive is removed after extraction when requested."""
        (path, extract, delete_kw) = archive
        extract(path, tmp_path / "out", **{delete_kw: True})
        assert not path.exists()


class TestTarValidity:
    """Tests for is_valid_tar and untar's handling of invalid archives."""

    @pytest.mark.parametrize("mode", TAR_MODES.values(), ids=TAR_MODES.keys())
    def test_valid(self, mode, tmp_path):
        """A well-formed archive of any compression is valid."""
        path = _make_tar(tmp_path, mode)
        assert is_valid_tar(path)
        assert is_valid_tar(str(path))

    def test_invalid(self, tmp_path):
        """A corrupt archive is invalid and untar refuses it."""
        bad = tmp_path / "bad.txz"
        bad.write_bytes(b"this is not a tar file")
        assert not is_valid_tar(bad)
        with pytest.raises(OSError):
            untar(bad, tmp_path / "out")
