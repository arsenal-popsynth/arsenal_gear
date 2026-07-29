"""
utils
=====
This subpackage contains various general purpose helper functions used across
the arsenal_gear package.
"""

from . import array_utils
from .citations import cite, find_citations, gather_bibtex
from .file_io import extract_one, find_match, get_metstr, is_valid_txz
from .math_utils import masked_power
from .scraper import downloader

__all__ = [
    "downloader",
    "is_valid_txz",
    "extract_one",
    "find_match",
    "masked_power",
    "array_utils",
    "get_metstr",
    "cite",
    "find_citations",
    "gather_bibtex",
]
