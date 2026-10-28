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
from io import BytesIO

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


def get_full_src_tree() -> bytes:
    """
    This function generates a tarball of the entire source tree returned by
    arsenal_gear.__path__ and returns it as a bytes object.  Subdirectories
    named 'data' or '__pycache__' are excluded.
    """
    buf = BytesIO()
    # we are writing the tarfile to memory in buf
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        tar.add(
            arsenal_gear.__path__[0],
            arcname="arsenal_gear",
            filter=_exclude_ignored_dirs,
        )
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


def construct_meta_dataset(f):
    """
    Add a metadata dataset to the h5py file object f.
    """
    ds = f.create_dataset("metadata/source", data=np.void(get_full_src_tree()))
    ds.attrs["format"] = "tar.gz"
    try:
        ds.attrs["git-commit"] = git_commit(arsenal_gear.__path__)
        ds.attrs["git-dirty"] = git_is_dirty(arsenal_gear.__path__)
    except subprocess.CalledProcessError:
        ds.attrs["git-commit"] = ""
