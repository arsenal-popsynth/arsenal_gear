"""
utils
==========

This file conatins various general purpose helper functions for the stellar evolution
module.
"""

import numpy as np

def _make_monotonic(x: np.float64, y:np.float64) -> np.float64:
    """
    Make an array montonically increasing by replacing decreasing values with
    interpolated values based on the nearest by points j and k such that 
    x[j] < x[k] and y[j] < y[k].
    Args:
        x (np.float64): indepednent variable, used for interpolation
        y (np.float64): dependent variable to be made monotonic
    """
    j = 0
    while(j < len(y)-1):
        if y[j] > y[j+1]:
            k = j+1
            while(y[j]>y[k]):
                k+=1
            y[j+1:k] = np.interp(x[j+1:k],[x[j],x[k]],[y[j],y[k]])
            j = k
        else:
            j+=1

def _index_monotonic(y: np.float64) -> np.ndarray:
    """
    Return the indices of the array y such that y is monotonic increasing
    """
    j = 0
    out = []
    while(j < len(y)-1):
        out.append(j)
        if y[j] > y[j+1]:
            k = j+1
            while(y[j]>y[k]):
                k += 1
            j = k+1
        else:
            j += 1
    return np.array(out, dtype=int)