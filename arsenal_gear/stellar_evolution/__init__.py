"""
stellar_evolution
=================

Submodule to handle stellar evolution codes and their isochrones
"""

from abc import ABC, abstractmethod

from . import isochrone

__all__ = ["isochrone"]
