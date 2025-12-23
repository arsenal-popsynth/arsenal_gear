"""
math_utils.py
=============

This file contains various utility functions for mathematical operations.
Functions:
    _masked_power: Computes power while handling masked values appropriately.
"""

from astropy.utils.masked import Masked
import numpy as np


def masked_power(base, exponent):
    """
    Computes base raised to the exponent, handling masked values appropriately.
    Args:
        base: The base value (can be masked).
        exponent: The exponent value (can be masked).
    Returns:
        The result of base ** exponent, with masked values preserved as masked.
    """
    if np.isscalar(base) and not np.isscalar(exponent):
        mask = Masked(exponent).mask
        nmask = np.logical_not(mask)
        res = np.zeros_like(exponent)
        res[nmask] = np.power(base, exponent[nmask])
        res = Masked(res,mask=mask)
    elif not(np.isscalar(base)) and np.isscalar(exponent):
        mask = Masked(base).mask
        nmask = np.logical_not(mask)
        res = np.zeros_like(base)
        res[nmask] = np.power(base[nmask], exponent)
        res = Masked(res,mask=mask)
    elif not np.isscalar(base) and not np.isscalar(exponent):
        if base.shape != exponent.shape:
            raise ValueError("Base and exponent must have the same shape.")
        mask = np.logical_or(Masked(base).mask, Masked(exponent).mask)
        nmask = np.logical_not(mask)
        res = np.zeros_like(base)
        res[nmask] = np.power(base[nmask], exponent[nmask])
        res = Masked(res,mask=mask)
    else:
        raise ValueError("Masking not needed for scalar power.")
    return res