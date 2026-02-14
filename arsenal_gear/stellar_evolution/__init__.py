"""
stellar_evolution
=================

Submodule to handle stellar evolution codes and their isochrones
"""

from . import isochrone
from ..utils.file_io import get_metstr

__all__ = ["Evolution"]

class Evolution():
    def __init__(self, **kwargs):
        # initialize the single star evolution interface
        se_opt = kwargs.get("single_evolution", "interpolator")
        self.mets = kwargs["metallicities"]
        self.nmets = len(self.mets)
        opt_list = isinstance(se_opt, list)
        opt_str = isinstance(se_opt, str)
        opt_iso = isinstance(se_opt, isochrone.AbstractIsochrone)
        self.ses = []
        if opt_list:
            # proces list of options
            if len(se_opt) != self.nmets:
                err_msg = "Length of single_evolution list must match \
                           number of metallicities."
                raise ValueError(err_msg)
            else:
                for i in range(self.nmets):
                    if isinstance(se_opt[i], isochrone.AbstractIsochrone):
                        if se_opt[i].met != self.mets[i]:
                            err_msg = f"Metallicity of single_evolution isochrone \
                                       at index {i} does not match metallicity \
                                       at index {i} in metallicities."
                            raise ValueError(err_msg)
                        else:
                            self.ses.append(se_opt[i])
                    elif isinstance(se_opt[i], str):
                        self.ses.append(self.handle_iso_str_single(se_opt[i], **kwargs))
                    else:
                        err_msg = f"Invalid single_evolution option type at index {i}: \
                                   {type(se_opt[i])}"
                        raise ValueError(err_msg)
        elif opt_str:
            if self.nmets > 1:
                self.ses = []
                for met in self.mets:
                    kwargs["met"] = met
                    self.ses.append(self.handle_iso_str_single(se_opt, **kwargs))
            else:
                kwargs["met"] = self.mets[0]
                self.ses.append(self.handle_iso_str_single(se_opt, **kwargs))
        elif opt_iso:
            self.ses.append(se_opt)
        else:
            err_msg = f"Invaliqd single_evolution option type: {type(se_opt)}"
            raise ValueError(err_msg)

        # initialize the binary star evolution interface
        self.be = None


    @staticmethod
    def handle_iso_str_single(opt:str, **kwargs) -> isochrone.AbstractIsochrone:
        if opt == "interpolator":
            return isochrone.IsochroneInterpolator(**kwargs)
        else:
            err_msg = f"Invalid single_evolution option string: {opt}"
            raise ValueError(err_msg)