"""
metadata.py
===========

Generate run-specific metadata to be written out into any arsenal_gear hdf5
binary outputs.
"""

from __future__ import annotations

import importlib
import importlib.util
import pkgutil
import tarfile
from io import BytesIO

import arsenal_gear


def subpackage_paths(pkg_name):
    """Map each package under pkg_name (including itself) to its __path__ entries."""
    pkg = importlib.import_module(pkg_name)
    result = {pkg_name: list(pkg.__path__)}

    for info in pkgutil.walk_packages(pkg.__path__, prefix=pkg_name + "."):
        if info.ispkg:
            spec = importlib.util.find_spec(info.name)
            result[info.name] = list(spec.submodule_search_locations)
    return result


def get_full_src_tree() -> bytes:
    """
    This function generates a tarball of the entire source tree returned by
    arsenal_gear.__path__ and returns it as a bytes object.
    """
    buf = BytesIO()
    # we are writing the tarfile to memory in buf
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        tar.add(arsenal_gear.__path__, arcname="dir")
    return buf.getvalue()
