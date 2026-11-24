"""
extract_source.py
=================

Recover the arsenal_gear source tree that was embedded in an HDF5 output file
by :func:`arsenal_gear.file_io.metadata.construct_meta_dataset`.

Usage::

    arsenal-extract-source output.h5              # extracts to ./output_source/
    arsenal-extract-source output.h5 -o src/      # extracts to ./src/
    arsenal-extract-source output.h5 --info       # print provenance, don't extract
"""

from __future__ import annotations

import argparse
import io
import os
import sys
import tarfile

import h5py

SOURCE_DATASET = "metadata/source"


def read_provenance(f: h5py.File) -> dict:
    """
    Collect the run provenance stored in the metadata group of an open
    arsenal_gear output file.  Fields that are absent are left out.
    """
    info = {}
    for key in ("rundate", "hostname"):
        if f"metadata/{key}" in f:
            info[key] = f[f"metadata/{key}"].asstr()[()]
    if SOURCE_DATASET in f:
        attrs = f[SOURCE_DATASET].attrs
        info["git-commit"] = attrs.get("git-commit", "") or "unknown"
        info["git-dirty"] = bool(attrs.get("git-dirty", True))
    return info


def read_source_blob(f: h5py.File) -> bytes:
    """
    Return the raw tar.gz bytes of the embedded source tree.
    """
    if SOURCE_DATASET not in f:
        raise KeyError(f"no '{SOURCE_DATASET}' dataset in {f.filename}")
    ds = f[SOURCE_DATASET]
    fmt = ds.attrs.get("format", "tar.gz")
    if fmt != "tar.gz":
        raise ValueError(f"unsupported source archive format '{fmt}'")
    return ds[()].tobytes()


def _safe_members(tar: tarfile.TarFile, dest: str):
    """
    Fallback for Pythons without tarfile extraction filters: refuse any member
    that would land outside dest or isn't a plain file or directory.
    """
    dest = os.path.realpath(dest)
    for member in tar.getmembers():
        target = os.path.realpath(os.path.join(dest, member.name))
        if os.path.commonpath([dest, target]) != dest:
            raise ValueError(f"refusing to extract {member.name!r} outside {dest}")
        if not (member.isfile() or member.isdir()):
            raise ValueError(f"refusing to extract non-regular member {member.name!r}")
        yield member


def extract_source(blob: bytes, dest: str) -> list[str]:
    """
    Extract the tar.gz bytes in blob into the directory dest and return the
    names of the extracted members.
    """
    os.makedirs(dest, exist_ok=True)
    with tarfile.open(fileobj=io.BytesIO(blob), mode="r:gz") as tar:
        if hasattr(tarfile, "data_filter"):
            tar.extractall(dest, filter="data")
        else:  # Python < 3.10.12 / 3.11.4
            tar.extractall(dest, members=_safe_members(tar, dest))
        return tar.getnames()


def main(argv: list[str] | None = None) -> int:
    """
    Command line entry point for arsenal-extract-source.
    """
    parser = argparse.ArgumentParser(
        description="Recover the arsenal_gear source tree embedded in an HDF5 output."
    )
    parser.add_argument("h5file", help="arsenal_gear HDF5 output file")
    parser.add_argument(
        "-o",
        "--output",
        help="directory to extract into (default: <h5file stem>_source)",
    )
    parser.add_argument(
        "--info",
        action="store_true",
        help="print the run provenance and exit without extracting",
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="extract even if the output directory already has files in it",
    )
    args = parser.parse_args(argv)

    try:
        with h5py.File(args.h5file, "r") as f:
            info = read_provenance(f)
            blob = None if args.info else read_source_blob(f)
    except (OSError, KeyError, ValueError) as err:
        print(f"error: {err}", file=sys.stderr)
        return 1

    for key, value in info.items():
        print(f"{key:>10}: {value}")
    if info.get("git-dirty"):
        print(
            "   warning: the source had uncommitted changes (or no git commit "
            "could be found), so the commit hash alone won't reproduce it."
        )
    if args.info:
        return 0

    dest = args.output or os.path.splitext(os.path.basename(args.h5file))[0] + "_source"
    if os.path.isdir(dest) and os.listdir(dest) and not args.force:
        print(
            f"error: {dest} already exists and is not empty (use --force to "
            "extract into it anyway)",
            file=sys.stderr,
        )
        return 1

    try:
        names = extract_source(blob, dest)
    except (tarfile.TarError, ValueError) as err:
        print(f"error: could not extract source: {err}", file=sys.stderr)
        return 1
    print(f"Extracted {len(names)} entries to {dest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
