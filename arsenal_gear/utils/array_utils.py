"""
array_utils.py
==============

This file contains various utility manipulation functions for arrays.
Functions:
    downloader: Method for downloading data from the web
"""

import numpy as np


def make_monotonic_increasing(x: np.float64, y: np.float64) -> np.float64:
    """
    Make an array montonically increasing by replacing decreasing values with
    interpolated values based on the nearest by points j and k such that
    x[j] < x[k] and y[j] < y[k].
    Args:
        x (np.float64): indepednent variable, used for interpolation
        y (np.float64): dependent variable to be made monotonic
    """
    j = 0
    n = len(y)
    while j < n-1:
        if y[j] >= y[j+1]:
            k = j+1
            while((k < n) and (y[j]>=y[k])):
                k+=1
            if k == n:
                k = n - 1
                if y[j] > y[k]:
                    last = y[j]*1.0001
                else:
                    last = y[k]
                y[j + 1 :] = np.interp(x[j + 1 :], [x[j], x[k]], [y[j], last])
            else:
                y[j + 1 : k] = np.interp(x[j + 1 : k], [x[j], x[k]], [y[j], y[k]])
            j = k
        else:
            j += 1
    return y


def make_monotonic_decreasing(x: np.float64, y: np.float64) -> np.float64:
    """
    Make an array montonically decreasing by replacing increasing values with
    interpolated values based on the nearest by points j and k such that
    x[j] < x[k] and y[j] > y[k].
    Args:
        x (np.float64): indepednent variable, used for interpolation
        y (np.float64): dependent variable to be made monotonic
    """
    j = 0
    n = len(y)
    while j < n - 1:
        if y[j] < y[j + 1]:
            k = j + 1
            while (k < n) and (y[j] < y[k]):
                k += 1
            if k == n:
                k = n - 1
                if y[j] < y[k]:
                    last = y[j] * 0.9
                else:
                    last = y[k]
                y[j + 1 :] = np.interp(x[j + 1 :], [x[j], x[k]], [y[j], last])
            else:
                y[j + 1 : k] = np.interp(x[j + 1 : k], [x[j], x[k]], [y[j], y[k]])
            j = k
        else:
            j += 1
    return y


def index_monotonic(y: np.float64) -> np.ndarray:
    """
    Return the indices inds, of the array y such that y[inds] is monotonic increasing
    """
    j = 0
    out = []
    while j < len(y) - 1:
        out.append(j)
        if y[j] > y[j + 1]:
            k = j + 1
            while y[j] > y[k]:
                k += 1
            # this should technically be k not k+1, but that
            # leads to weird behavior, possibly due to very small
            # differences in values, so skipping to the next index seems safer
            j = k + 1
        else:
            j += 1
    return np.array(out, dtype=int)
