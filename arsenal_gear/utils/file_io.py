"""
file_io.py
==========

This file contains various utility functions for file input/output operations,
including downloading and extracting data from online sources.
All functions taking paths accept either strings or pathlib.Path objects.
Functions:
    downloader: Method for downloading data from the web
    untar: Extract (part of) a tar archive (.tar, .tar.gz, .txz, ...)
    unzip: Extract (part of) a zip archive
    is_valid_tar: Check whether a file is a readable tar archive
    find_match: Find files matching glob patterns
    get_metstr: Convert a metallicity to a MIST-style string
"""

import glob
import os.path as osp
import tarfile
from pathlib import Path
from zipfile import ZipFile

import requests
from tqdm import tqdm


def downloader(fname, url, message=None, username=None, password=None):
    """
    Method for downloading data from the web where available.

    Args:
        fname (str or Path): The file name or path to save the downloaded file.
        url (str): The URL of the file to download.
        message (str): Optional message to display before downloading.
        username (str): Optional username for password-protected downloads.
        password (str): Optional password for password-protected downloads.

    Raises:
        Exception: If the download fails.

    """
    if message is not None:
        print(message)

    if username is not None:
        auth_info = (username, password)
    else:
        auth_info = None

    try:
        response = requests.get(
            url,
            stream=True,
            timeout=10,
            auth=auth_info,
        )
        response.raise_for_status()
    except requests.exceptions.Timeout as e:
        raise TimeoutError("Request timed out. Check internet connection.") from e
    except requests.exceptions.ConnectionError as e:
        raise ConnectionError("Connection error. Check internet connection.") from e
    except requests.exceptions.HTTPError as e:
        raise RuntimeError(f"HTTP error occurred: {e}") from e
    except requests.exceptions.TooManyRedirects as e:
        raise RuntimeError("Too many redirects. Check the URL.") from e
    except requests.exceptions.RequestException as e:
        raise RuntimeError(f"Download failed: {e}") from e

    # Get file size
    total_size = int(response.headers.get("content-length", 0))
    # create a progress bar
    tqdm_args = {
        "desc": "Downloading",
        "total": total_size,
        "unit": "B",
        "unit_scale": True,
        "unit_divisor": 1024,
    }
    # write the file
    with open(fname, "wb") as f, tqdm(**tqdm_args) as prog_bar:
        for chunk in response.iter_content(chunk_size=1024):
            f.write(chunk)
            prog_bar.update(len(chunk))


def untar(
    tar_file, target_dir=None, target_file=None, delete_tar=False, inspect=False
) -> Path:
    """
    Un-compress a tar file. Compression (none, gz, xz, bz2) is detected automatically.

    Args:
        tar_file    (str or Path): Path to the tar file
        target_dir  (str or Path): Directory in which to un-compress the tar file
                                   If None, use the directory containing the tar file
        target_file (str): Name of member to extract from the tar file
                           If None, extract all members
        delete_tar (bool): Delete the tar file after extracting
        inspect    (bool): Print the names of the members in the tar file
                           without extracting anything
    Returns:
        target_dir (Path): Directory into which the tar file was extracted
    """
    fname = Path(tar_file)
    target_dir = fname.parent if target_dir is None else Path(target_dir)

    if not tarfile.is_tarfile(fname):
        raise OSError(
            f"{fname} is not a valid tar file. Try again with `force_download=True`"
        )
    with tarfile.open(fname, "r:*") as tar:
        if inspect:
            for name in tar.getnames():
                print(name)
            return target_dir
        target_dir.mkdir(parents=True, exist_ok=True)
        print(f"Extracting {fname}...")
        if target_file is None:
            tar.extractall(path=target_dir)
        elif target_file in tar.getnames():
            tar.extract(target_file, path=target_dir)
        else:
            raise FileNotFoundError(f"{target_file} not found in {fname}")

    if delete_tar and fname.exists():
        print("Deleting", fname)
        fname.unlink()

    return target_dir


def unzip(
    zip_file, target_dir=None, target_file=None, delete_zip=False, inspect=False
) -> Path:
    """
    Un-compress a zip file.

    Args:
        zip_file    (str or Path): Path to the zip file
        target_dir  (str or Path): Directory in which to un-compress the zip file
                                   If None, use the directory containing the zip file
        target_file (str): Name of member to extract from the zip file
                           If None, extract all members
        delete_zip (bool): Delete the zip file after extracting
        inspect    (bool): Print the names of the members in the zip file
                           without extracting anything
    Returns:
        target_dir (Path): Directory into which the zip file was extracted
    """
    fname = Path(zip_file)
    target_dir = fname.parent if target_dir is None else Path(target_dir)

    with ZipFile(fname, "r") as zip_archive:
        file_names = zip_archive.namelist()
        if inspect:
            for name in file_names:
                print(name)
            return target_dir
        target_dir.mkdir(parents=True, exist_ok=True)
        if target_file is None:
            zip_archive.extractall(path=target_dir)
        elif target_file in file_names:
            zip_archive.extract(target_file, path=target_dir)
        else:
            raise FileNotFoundError(f"{target_file} not found in {fname}")

    if delete_zip and fname.exists():
        print("Deleting", fname)
        fname.unlink()

    return target_dir


def is_valid_tar(fname) -> bool:
    """
    Check if a tar file (of any compression type) is valid.

    Parameters:
    - fname (str or Path): The path to the tar file.

    Returns:
    - bool: True if the tar file is valid, False otherwise.
    """
    try:
        with tarfile.open(fname, "r:*") as tar:
            tar.getmembers()
        return True
    except (tarfile.TarError, ValueError, OSError, EOFError) as e:
        print(f"Invalid tar file: {e}")
        return False


def get_metstr(met: float) -> str:
    """
    Convert a metallicity value to a string format suitable for file naming.

    Parameters:
    - met (float): The metallicity value to convert, assumed to be in log10(Z/Zsun)
                   (or at least on a log scale)

    Returns:
    - str: A string representation of the metallicity
           'p' indicates a postive metallicity (super-solar)
           'm' indicates a negative metallicity (sub-solar)
    """
    if met < 0:
        metstr = f"m{-1*met:.2f}"
    else:
        metstr = f"p{met:.2f}"
    return metstr


def find_match(patterns, base_dir=None):
    """
    Finds files in base_dir matching one of the glob patterns in patterns.
    """
    if base_dir is None:
        glob_match = lambda p: sorted(glob.glob(*p))
    else:
        glob_match = lambda p: sorted(glob.glob(osp.join(base_dir, *p)))
    for p in patterns:
        f = glob_match(p)
        if f:
            break
    return f
