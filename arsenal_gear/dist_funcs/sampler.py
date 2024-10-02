"""
sampler
=======

This submodule contains the base class for anything that will draw a sample from a ProbDistFunc.
"""

from typing import Type

from . import ProbDistFunc


class Sampler():
    """
    This class is the superclass of all samplers (IMF, Binary Mass Function,
    Binary Radius function, etc.) that arsenal will use.

    :param PDF: the probability distribution function we will draw from
    :type PDF: ProbDistFunc
    """

    def __init__(self, PDF:Type(ProbDistFunc)) -> None:
        self.PDF = PDF

    def __call__(self, N:int) -> np.float64:
        """
        Return an N-element numpy array of elements drawn from the PDF.

        :param N: the number of elements drawn
        :type N: int
        :rtype: np.float64
        """
        pass
