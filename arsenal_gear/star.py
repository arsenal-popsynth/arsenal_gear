"""
star
====

This submodule defines the class for individual stars.
"""

from typing import Type

__all__ = ['star']

class Star():
    def __init__(self, mass: float, metals: float, age: float=0) -> None:
        self.mass = mass
        self.metals = metals
        self.age = age

