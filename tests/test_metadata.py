"""
test_metadata
==========

This file contains unit tests for the run metadata written into hdf5 outputs,
which lives in arsenal_gear/file_io/metadata.py
"""

# A test that takes a fixture names its parameter after that fixture, which
# pylint cannot tell apart from shadowing a global.
# pylint: disable=redefined-outer-name

import io
import shutil
import subprocess
import tarfile
from datetime import datetime
from socket import gethostname

import h5py
import pytest

from arsenal_gear.file_io import metadata
from arsenal_gear.file_io.metadata import (
    IGNORED_DIRS,
    _exclude_ignored_dirs,
    construct_meta_dataset,
    get_full_src_tree,
    git_commit,
    git_is_dirty,
)

needs_git = pytest.mark.skipif(shutil.which("git") is None, reason="git not found")


@pytest.fixture
def h5file(tmp_path):
    """A fresh, writable hdf5 file."""
    with h5py.File(tmp_path / "test.h5", "w") as f:
        yield f


def _git(repo, *args):
    """Run a git command in repo with a throwaway identity."""
    subprocess.run(
        ["git", "-c", "user.name=test", "-c", "user.email=test@test", *args],
        cwd=repo,
        check=True,
        capture_output=True,
    )


@pytest.fixture
def git_repo(tmp_path):
    """A git repository with a single clean commit."""
    repo = tmp_path / "repo"
    repo.mkdir()
    _git(repo, "init")
    (repo / "file.txt").write_text("hello\n")
    _git(repo, "add", "file.txt")
    _git(repo, "commit", "-m", "initial")
    return repo


# --- the tarball filter ---------------------------------------------------


def _tarinfo(name, is_dir):
    info = tarfile.TarInfo(name)
    info.type = tarfile.DIRTYPE if is_dir else tarfile.REGTYPE
    return info


@pytest.mark.parametrize("dirname", IGNORED_DIRS)
def test_filter_drops_ignored_dirs(dirname):
    assert _exclude_ignored_dirs(_tarinfo(f"arsenal_gear/sub/{dirname}", True)) is None


def test_filter_keeps_other_dirs():
    info = _tarinfo("arsenal_gear/file_io", True)
    assert _exclude_ignored_dirs(info) is info


@pytest.mark.parametrize("name", ["arsenal_gear/data", "arsenal_gear/data.py"])
def test_filter_keeps_files_named_like_ignored_dirs(name):
    """Only directories are filtered, not files that share their name."""
    info = _tarinfo(name, False)
    assert _exclude_ignored_dirs(info) is info


def test_filter_matches_whole_name_only():
    info = _tarinfo("arsenal_gear/metadata", True)
    assert _exclude_ignored_dirs(info) is info


# --- the source tarball ---------------------------------------------------


@pytest.fixture(scope="module")
def src_members():
    """Member names of the source tarball."""
    with tarfile.open(fileobj=io.BytesIO(get_full_src_tree()), mode="r:gz") as tar:
        return tar.getnames()


def test_src_tree_is_gzipped_bytes():
    blob = get_full_src_tree()
    assert isinstance(blob, bytes)
    assert blob[:2] == b"\x1f\x8b"  # gzip magic number


def test_src_tree_contains_source(src_members):
    assert "arsenal_gear/__init__.py" in src_members
    assert "arsenal_gear/file_io/metadata.py" in src_members


def test_src_tree_is_rooted_at_package(src_members):
    assert all(
        name == "arsenal_gear" or name.startswith("arsenal_gear/")
        for name in src_members
    )


def test_src_tree_excludes_ignored_dirs(src_members):
    for name in src_members:
        assert not set(name.split("/")) & set(IGNORED_DIRS), name


# --- git state ------------------------------------------------------------


@needs_git
def test_git_commit_returns_head(git_repo):
    head = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=git_repo,
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    commit = git_commit(git_repo)
    assert commit == head
    assert len(commit) == 40
    int(commit, 16)  # it's a hex sha


@needs_git
def test_git_clean_repo_is_not_dirty(git_repo):
    assert git_is_dirty(git_repo) is False


@needs_git
def test_git_modified_file_is_dirty(git_repo):
    (git_repo / "file.txt").write_text("changed\n")
    assert git_is_dirty(git_repo) is True


@needs_git
def test_git_untracked_file_is_dirty(git_repo):
    (git_repo / "new.txt").write_text("new\n")
    assert git_is_dirty(git_repo) is True


@needs_git
@pytest.mark.parametrize("func", [git_commit, git_is_dirty])
def test_git_outside_repo_raises(tmp_path, monkeypatch, func):
    # Stop git from walking up into whatever repo tmp_path happens to live in.
    monkeypatch.setenv("GIT_CEILING_DIRECTORIES", str(tmp_path))
    with pytest.raises(subprocess.CalledProcessError):
        func(tmp_path)


# --- the metadata group ---------------------------------------------------


def test_meta_dataset_contents(h5file):
    before = datetime.now()
    construct_meta_dataset(h5file, "paramfile contents")
    after = datetime.now()

    meta = h5file["metadata"]
    rundate = datetime.fromisoformat(meta["rundate"].asstr()[()])
    assert before <= rundate <= after
    assert meta["hostname"].asstr()[()] == gethostname()
    assert meta["paramfile"].asstr()[()] == "paramfile contents"
    assert "bibtex" not in meta


def test_meta_dataset_bibtex(h5file):
    bibtex = "@ARTICLE{key,\n  title = {Entry}\n}"
    construct_meta_dataset(h5file, "", bibtex=bibtex)
    assert h5file["metadata/bibtex"].asstr()[()] == bibtex


def test_meta_dataset_source_roundtrips(h5file):
    construct_meta_dataset(h5file, "")
    ds = h5file["metadata/source"]
    assert ds.attrs["format"] == "tar.gz"
    blob = ds[()].tobytes()
    with tarfile.open(fileobj=io.BytesIO(blob), mode="r:gz") as tar:
        assert "arsenal_gear/file_io/metadata.py" in tar.getnames()


def test_meta_dataset_git_attrs(h5file, monkeypatch):
    monkeypatch.setattr(metadata, "git_commit", lambda path: "abc123")
    monkeypatch.setattr(metadata, "git_is_dirty", lambda path: True)
    construct_meta_dataset(h5file, "")
    attrs = h5file["metadata/source"].attrs
    assert attrs["git-commit"] == "abc123"
    assert attrs["git-dirty"]


@pytest.mark.parametrize(
    "error",
    [
        subprocess.CalledProcessError(128, ["git", "rev-parse", "HEAD"]),
        FileNotFoundError("git"),
    ],
    ids=["not-a-repo", "no-git"],
)
def test_meta_dataset_without_git(h5file, monkeypatch, error):
    """An install outside a git checkout, or without git, still writes
    metadata, and doesn't claim its source is clean."""

    def fail(path):
        raise error

    monkeypatch.setattr(metadata, "git_commit", fail)
    construct_meta_dataset(h5file, "")
    attrs = h5file["metadata/source"].attrs
    assert attrs["git-commit"] == ""
    assert attrs["git-dirty"]


@needs_git
def test_meta_dataset_real_git(h5file):
    """Use the real git helpers against the arsenal_gear source tree, without
    mocking, to check they're called with a usable path."""
    construct_meta_dataset(h5file, "")
    assert "git-commit" in h5file["metadata/source"].attrs
