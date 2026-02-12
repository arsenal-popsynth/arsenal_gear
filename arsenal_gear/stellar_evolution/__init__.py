"""
stellar_evolution
=================

Submodule to handle stellar evolution codes and their isochrones
"""

from . import isochrone

__all__ = ["Evolution"]

class Evolution():
    def __init__(self, **kwargs):

        # initialize the single star evolution interface
        self.se = kwargs.get("single_evolution", "interpolator")
        if not isinstance(self.se, isochrone.AbstractIsochrone):
            if self.se == "interpolator":
                self.se = isochrone.IsochroneInterpolator(**kwargs)
            else:
                err_msg = f"Invalid single_evolution option: {self.se}"
                raise ValueError(err_msg)
        
        # initialize the binary star evolution interface
        self.be = None
