import astropy.units as u
import numpy as np
import pytest
from numpy.testing import assert_array_equal

import arsenal_gear


def test_limits():
    pdf = arsenal_gear.dist_funcs.ProbDistFunc(0, 1)
    assert pdf(np.array([-1])) == 0
    assert pdf(np.array([2])) == 0
    assert pdf(np.array([0])) == 1
    assert pdf(np.array([1])) == 1
    assert pdf(np.array([0.5])) == 1

def test_normalization():
    pdf = arsenal_gear.dist_funcs.ProbDistFunc(0, 10, normalized=True)
    assert pdf(np.array([-1])) == 0
    assert pdf(np.array([2])) == 0.1
    assert pdf(np.array([10])) == 0.1
