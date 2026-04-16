"""
formation
=================

Submodule to describe the properties of a stellar population
when it is first formed.
"""

import astropy.units as u
import numpy as np

# local imports
from . import dist_funcs
from .form_data_structures import SinglePop, BinaryPop

__all__ = ["Formation", "SinglePop", "BinaryPop"]


class Formation():
    """
    Abstract base class for describing the properties of a stellar population
    at the time of formation.
    """
    def __init__(self, **kwargs) -> None:
        # get population parameters from kwargs
        subpop_dicts = [kwargs.get(k) for k in kwargs if k.startswith("pop")]
        self.nsubpops = len(subpop_dicts)
        if (self.nsubpops == 0):
            default_pop = {
                "Mtot": 1e6 * u.Msun,
                "metallicity": 0.0 * u.dimensionless_unscaled,
                "rotation": 0.0 * u.dimensionless_unscaled,
                "imf": "salpeter",
                "mmin": 0.08,
                "mmax": 100.0,
                "discrete": True,
                "seed": kwargs.get("seed", None)}
            subpop_dicts = [default_pop]
            self.nsubpops = 1
        for i, pop_dict in enumerate(subpop_dicts):
            if not isinstance(pop_dict, dict):
                err_msg = (f"Subpopulation {i} is not a dictionary. "
                           "Please provide subpopulation parameters as dictionaries."
                )
                raise ValueError(err_msg)
        
        subpops = []
        for pop_dict in subpop_dicts:
            pop_type = pop_dict.get("type", "single")
            if pop_type == "single":
                subpops.append(self._init_single_pop(pop_dict))
            elif pop_type == "binary":
                subpops.append(self._init_binary_pop(pop_dict))
            else:
                err_msg = (f"Invalid population type: {pop_type}. "
                           "Must be 'single' or 'binary'."
                )
                raise ValueError(err_msg)
        self.subpops = subpops
        # Nstar can be a float if one of the populations is not discrete
        self.tforms = np.array([pop.tform.to(u.Myr).value for pop in self.subpops])*u.Myr
        self.metallicities = [pop.metallicity for pop in self.subpops]


    @staticmethod
    def _init_single_pop(pop_dict: dict) -> SinglePop:
        """
        Initialize a single population of stars from a dictionary of parameters.
        """
        # check Mtot value
        Mtot = pop_dict.get("Mtot", 1e6 * u.Msun)
        if isinstance(Mtot, u.Quantity):
            Mtot = Mtot.to(u.Msun)
        elif isinstance(Mtot, (int, float)):
            Mtot = Mtot * u.Msun
        else:
            err_msg = (f"Invalid Mtot: {Mtot}. Must be a quantity specifying "
                       "mass or a number specifying mass in Msun."
            )
            raise ValueError(err_msg)
        
        # check tform value
        tform = pop_dict.get("tform", 0.0 * u.Myr)
        if isinstance(tform, u.Quantity):
            tform = tform.to(u.Myr)
        elif isinstance(tform, (int, float)):
            tform = tform * u.Myr
        else:
            err_msg = (f"Invalid tform: {tform}. Must be a quantity specifying "
                       "time or a number specifying time in Myr."
            )
            raise ValueError(err_msg)

        # check metallicity
        metallicity = pop_dict.get("metallicity", 0.0 * u.dimensionless_unscaled)
        if isinstance(metallicity, u.Quantity):
            metallicity = metallicity.to(u.dimensionless_unscaled)
        elif isinstance(metallicity, (int, float)):
            metallicity = metallicity * u.dimensionless_unscaled
        else:
            err_msg = (f"Invalid metallicity: {metallicity}. Must be a dimensionless quantity "
                       "or a number specifying metallicity in log10(Z/Zsun)."
            )
            raise ValueError(err_msg)
        
        # check rotation
        rotation = pop_dict.get("rotation", 0.0 * u.dimensionless_unscaled)
        if isinstance(rotation, u.Quantity):
            rotation = rotation.to(u.dimensionless_unscaled)
        elif isinstance(rotation, (int, float)):
            rotation = rotation * u.dimensionless_unscaled
        else:
            err_msg = (f"Invalid rotation: {rotation}. Must be a dimensionless quantity "
                       "or a number specifying rotation as a fraction of critical rotation."
            )
            raise ValueError(err_msg)
        
        # check IMF
        imf = pop_dict.get("imf", "salpeter")
        if isinstance(imf, dist_funcs.imf.IMF):
            mmin = imf.min_mass * u.Msun
            mmax = imf.max_mass * u.Msun
        elif isinstance(imf, str):
            mmin = pop_dict.get("mmin", 0.08)
            mmax = pop_dict.get("mmax", 100.)
            try:
                if isinstance(mmin, u.Quantity):
                    mmin = mmin.to(u.Msun).value
                if isinstance(mmax, u.Quantity):
                    mmax = mmax.to(u.Msun).value
                (mmin, mmax) = (float(mmin), float(mmax))
            except Exception as exc:
                err_msg = (f"Invalid mmin or mmax: {mmin}, {mmax}. "
                           "Must be floats or quantities specifying mass in Msun."
                )
                raise ValueError(err_msg) from exc
            mmin *= u.Msun
            mmax *= u.Msun
            seed = pop_dict.get("seed", None)
            if imf == "salpeter":
                imf = dist_funcs.imf.Salpeter(mmin, mmax, seed=seed)
            elif imf == "kroupa":
                imf = dist_funcs.imf.Kroupa(mmin, mmax, seed=seed)
            elif imf == "miller-scalo":
                imf = dist_funcs.imf.MillerScalo(mmin, mmax, seed=seed)
            elif imf == "chabrier":
                imf = dist_funcs.imf.Chabrier(mmin, mmax, seed=seed)
            else:
                err_msg = (f"Invalid IMF string specified, {imf}."
                            " Currently supported options are 'salpeter', 'kroupa',"
                            "'miller-scalo', and 'chabrier'."
                )
                raise ValueError(err_msg)
        else:
            err_msg = (f"Invalid IMF specified, {imf}. "
                       "Must be an IMF instance or string specifying IMF."
            )
            raise ValueError(err_msg)

        Nstar = (Mtot / imf.mean()).value

        discrete = pop_dict.get("discrete", True)
        if not isinstance(discrete, bool):
            err_msg = f"Invalid discrete parameter: {discrete}. Must be a boolean."
            raise ValueError(err_msg)
        if discrete:
            masses = imf.sample_mass(Mtot)
            Nstar = len(masses)
        else:
            masses = None

        return SinglePop(Mtot,
                         tform,
                         Nstar,
                         metallicity,
                         rotation,
                         imf,
                         mmin,
                         mmax,
                         discrete,
                         masses)
    
    @staticmethod
    def _init_binary_pop(pop_dict: dict) -> BinaryPop:
        """
        Initialize a binary population of stars from a dictionary of parameters.
        """
        raise NotImplementedError("Binary populations not yet implemented.")

    @property
    def Mtot(self) -> u.Quantity["mass"]:
        """
        Return the total mass of the population, which is just the sum of Mtot
        for all subpopulations.
        """
        return np.sum([pop.Mtot.to(u.Msun).value for pop in self.subpops]) * u.Msun

    @property
    def Nstar(self) -> int | float:
        """
        Return the total number of stars in the population, which is just the
        sum of Nstar for all subpopulations.
        """
        Nstar = 0
        for pop in self.subpops:
            Nstar += pop.Nstar
        return Nstar

    @property
    def mean_mass(self) -> u.Quantity["mass"]:
        """
        Return the mean mass of the population, which is just Mtot/Nstar.
        """
        if self.Nstar == 0:
            return 0.0 * u.Msun
        else:
            return (self.Mtot / self.Nstar).to(u.Msun)
