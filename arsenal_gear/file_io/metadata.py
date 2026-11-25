"""
metadata.py
===========

Generate run-specific metadata to be written out into any arsenal_gear hdf5
binary outputs.
"""

from __future__ import annotations

import os
import subprocess
import tarfile
from dataclasses import dataclass
from datetime import datetime
from io import BytesIO
from socket import gethostname

import numpy as np

import arsenal_gear

IGNORED_DIRS = ("data", "__pycache__")


def _exclude_ignored_dirs(tarinfo: tarfile.TarInfo) -> tarfile.TarInfo | None:
    """
    Filter for tarfile.add that drops any subdirectory named in IGNORED_DIRS
    (along with everything beneath it).
    """
    if tarinfo.isdir() and os.path.basename(tarinfo.name) in IGNORED_DIRS:
        return None
    return tarinfo


def get_full_src_tree(path: str | os.PathLike | None = None) -> bytes:
    """
    This function generates a tarball of the entire source tree at path
    (by default, the directory in arsenal_gear.__path__) and returns it as a
    bytes object.  Subdirectories named 'data' or '__pycache__' are excluded.

    :param path: Directory to archive; defaults to the imported package
    :type path: str or os.PathLike, optional
    :return: Gzipped tarball of the source tree, with all paths rooted at
        ``arsenal_gear/``
    :rtype: bytes
    """
    if path is None:
        path = arsenal_gear.__path__[0]
    buf = BytesIO()
    # we are writing the tarfile to memory in buf
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        tar.add(path, arcname="arsenal_gear", filter=_exclude_ignored_dirs)
    return buf.getvalue()


def git_commit(path="."):
    return subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=path,
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()


def git_is_dirty(path="."):
    out = subprocess.run(
        ["git", "status", "--porcelain"],
        cwd=path,
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    return bool(out.strip())


@dataclass(frozen=True)
class SourceSnapshot:
    """
    The arsenal_gear source tree and its git state, captured together at one
    moment.

    :ivar tarball: Gzipped tarball from :func:`get_full_src_tree`
    :ivar git_commit: Commit hash of the source tree, or ``""`` if there is
        none (not a git checkout, or git isn't installed)
    :ivar git_dirty: True if the working tree had uncommitted changes, or if
        no commit could be found
    """

    tarball: bytes
    git_commit: str
    git_dirty: bool


def snapshot_source(path: str | os.PathLike) -> SourceSnapshot:
    """
    Capture the source tree at path along with its git commit and dirty
    state.

    :param path: Directory holding the arsenal_gear package
    :type path: str or os.PathLike
    :rtype: SourceSnapshot
    """
    tarball = get_full_src_tree(path)
    try:
        commit = git_commit(path)
        dirty = git_is_dirty(path)
    except (subprocess.CalledProcessError, FileNotFoundError):
        # Not a git checkout, or git isn't installed.  With no commit to
        # compare against, the source can't be vouched for as clean.
        commit, dirty = "", True
    return SourceSnapshot(tarball=tarball, git_commit=commit, git_dirty=dirty)


# The source as it was when arsenal_gear was imported.  arsenal_gear/__init__.py
# imports this module before any of its other submodules, so this matches the
# code that is running even if the files on disk are edited afterwards.
SOURCE_SNAPSHOT = snapshot_source(arsenal_gear.__path__[0])


def construct_meta_dataset(f, paramfile: str, bibtex=None) -> None:
    """
    Add a metadata group to the h5py file object f, recording when, where,
    and with what source and parameters the run was made.  This writes
    ``metadata/rundate``, ``metadata/hostname``, ``metadata/paramfile``,
    ``metadata/source`` (with ``format``, ``git-commit`` and ``git-dirty``
    attributes), and ``metadata/bibtex`` if one is given.  The source and its
    git state come from :data:`SOURCE_SNAPSHOT`, taken when arsenal_gear was
    imported, not from the files on disk when this is called.

    :param f: Open, writable hdf5 file to add the metadata to
    :type f: h5py.File
    :param paramfile: Full text of the parameter file used for the run
    :type paramfile: str
    :param bibtex: BibTeX for the models used in the run
    :type bibtex: str, optional
    """
    f.create_dataset("metadata/rundate", data=datetime.now().isoformat())
    f.create_dataset("metadata/hostname", data=gethostname())
    f.create_dataset("metadata/paramfile", data=paramfile)
    if bibtex is not None:
        f.create_dataset("metadata/bibtex", data=bibtex)

    snap = SOURCE_SNAPSHOT
    ds = f.create_dataset("metadata/source", data=np.void(snap.tarball))
    ds.attrs["format"] = "tar.gz"
    ds.attrs["git-commit"] = snap.git_commit
    ds.attrs["git-dirty"] = snap.git_dirty
