"""
file_io.py
==========

This file contains various utility functions for file input/output operations.
Functions:
    downloader: Method for downloading data from the web
"""
import glob
import os
import os.path as osp
import tarfile


def is_valid_txz(fname):
    """
    Check if a .txz file is valid.

    Parameters:
    - fname (str): The path to the .txz file.

    Returns:
    - bool: True if the .txz file is valid, False otherwise.
    """
    try:
        with tarfile.open(fname, "r:xz") as tar:
            tar.getmembers()
        return True
    except (tarfile.TarError, ValueError, OSError, EOFError) as e:
        print(f"Invalid .txz file: {e}")
        return False

def extract_one(fname, extractdir, delete_txz=False):
    """
    Unzips a single ZIP file.
    """
    # Ensure output directory exists
    os.makedirs(extractdir, exist_ok=True)
    if not tarfile.is_tarfile(fname):
        raise IOError(f'{fname} is not a valid txz file. '
                        'Try again with `force_download=True`')
    with tarfile.open(fname, 'r:xz') as tar:
        tar.extractall(path=extractdir)

    if delete_txz and fname.exists():
        print(fname)
        fname.unlink()

    return extractdir

def find_match(patterns, base_dir=None):
    """
    Finds files in base_dir matching one of the glob patterns in patterns.
    """
    if base_dir is None:
        glob_match = lambda p: sorted(glob.glob(*p))
    else:
        glob_match = lambda p: sorted(glob.glob(osp.join(base_dir,*p)))
    for p in patterns:
        f = glob_match(p)
        if f:
            break
    return f
