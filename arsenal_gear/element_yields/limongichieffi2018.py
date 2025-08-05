"""
limongichieffi2018
==========

This submodule contains all the code required to load yields from
Limongi & Chieffi (2018).
"""
import os
import re
from typing import List

import astropy.units as u
from astropy.units import Quantity
import numpy as np

from .yields import Yields, Source


class LimongiChieffi2018(Yields):
    """
    Include yields majority of elements (and isotopes) for massive stars 
    [13, 15, 20, 25, 30, 40, 60, 80, 120] Msun, with initial iron abundance 
    ([Fe/H]) [0, -1, -2, -3], and three different models for stellar 
    rotation [0, 150, 300] km/s.

    Download yields and place/link to directory LimongiChieffi2018 placed
    in the path of this file, e.g.,
    ln -s /path/to/downloaded/yields /path/to/arsenal_gear/element_yields/LimongiChieffi2018

    Yields available at http://orfeo.iaps.inaf.it/

    Reference: Limongi M & Chieffi A, 2018, ApJS, 237, 13L
    """

    def __init__(self, model: str = 'R') -> None:
        """
        Args:
            model: choise of model to load, see Limongi & Chieffi (2018) for details.

        Usage:
            >> lc2018 = arsenal_gear.element_yields.LimongiChieffi2018()
            >> mass = np.linspace(8, 120, 1000)*u.M_sun
            >> plt.plot(mass, yields.ccsn_yields('H', 
                                                 mass=mass, rot=0*u.km/u.s, 
                                                 metal = 1.345e-2*u.dimensionless_unscaled, 
                                                 interpolate='nearest'), 
                        '-', color=colors[-1])

        Attributes:
            model            Available models.
            mass             Tabulated masses.
            metal            Tabulated metallicities.
            feh              Tabulated [Fe/H].
            rot              Tabulated stellar rotation velocities.
            ccsn_max         Assumed minimum mass for direct collase to black hole.
            filedir          Directory of this file (used for relative path).
            yield_tablefile  Table filename (total yields).
            wind_tablefile   Table filename (wind yields).
            elements         Elements available in table.
            atomic_num       Atomic numbers of available elements.
            wind             Source object for stellar winds (massive stars)
            ccsn             Source object for core-collapse SNe
        """

        self.models = ['F', 'I', 'M', 'R']
        self.mass = np.array(
            [13.0, 15.0, 20.0, 25.0, 30.0, 40.0, 60.0, 80.0, 120.0])
        self.metal = np.array([3.236e-5, 3.2363e-4, 3.236e-3, 1.345e-2])
        self.feh = np.array([-3, -2, -1, 0])
        self.rot = np.array([0, 150, 300])
        self.ccsn_mmax = 25

        if model not in self.models:
            raise ValueError(f"Model {model} does not exist.")

        self.filedir = os.path.dirname(os.path.realpath(__file__))

        self.yield_tablefile = self.filedir + \
            f'/LimongiChieffi2018/tab_{model}/tab_yieldstot_ele_exp.dec'
        self.wind_tablefile = self.filedir + \
            f'/LimongiChieffi2018/tab_{model}/tabwind.dec'

        self.elements, self.atomic_num = self.get_element_list()

        # Stellar wind yields
        wind_yld = self.load_wind_yields()
        self.wind = Source(
            self.elements, [self.rot, self.metal, self.mass], wind_yld)

        # Core-collapse SNe yields
        ccsn_yld = self.load_ccsn_yields()
        self.ccsn = Source(
            self.elements, [self.rot, self.metal, self.mass], ccsn_yld)

    def ccsn_yields(self,
                    elements: List[str],
                    mass: Quantity["mass"],
                    metal: Quantity["dimensionless"],
                    rot: Quantity["velocity"],
                    interpolate: str = "nearest",
                    extrapolate: bool = False) -> Quantity["mass"]:
        """ Interpolate yields from core-collapse supernovae for specified elements. 
            Stellar parameters can be provided as single value, array + single value, or arrays.  

            Args:
                elements: list of elements, as specified by symbols (e.g., ['H'] for hydrogen).
                mass: stellar masses, single or list/array
                metal: stellar metallicity, single or list/array
                rot: stellar rotation, single or list/array
                interpolate: passed as method to scipy.interpolate.RegularGridInterpolator
                extrapolate: if False, then mass, metal, and rot are set to limits if outside bound. 
            Returns:
                List of yields matching provided element list

        """
        args = [rot.to(u.km/u.s).value, metal.value, mass.to(u.M_sun).value]
        return self.ccsn.get_yld(elements, args,
                                 interpolate=interpolate, extrapolate=extrapolate)*u.M_sun

    def wind_yields(self,
                    elements: List[str],
                    mass: Quantity["mass"],
                    metal: Quantity["dimensionless"],
                    rot: Quantity["velocity"],
                    interpolate: str = "nearest",
                    extrapolate: bool = False) -> Quantity["mass"]:
        """ Interpolate yields from massive stars ejected as winds for specified elements. 
            Stellar parameters can be provided as single value, array + single value, or arrays.  

            Args:
                elements: list of elements, as specified by symbols (e.g., ['H'] for hydrogen).
                mass: stellar masses, single or list/array
                metal: stellar metallicity, single or list/array
                rot: stellar rotation, single or list/array
                interpolate: passed as method to scipy.interpolate.RegularGridInterpolator
                extrapolate: if False, then mass, metal, and rot are set to table if outside bound. 
            Returns:
                List of yields matching provided element list

        """
        args = [rot.to(u.km/u.s).value, metal.value, mass.to(u.M_sun).value]
        return self.wind.get_yld(elements, args,
                                 interpolate=interpolate, extrapolate=extrapolate)*u.M_sun

    def get_element_list(self) -> None:
        """ Read element symbols and atomic numbers from tables.
        """
        with open(self.wind_tablefile, 'r', encoding="utf-8") as file:
            lines = file.readlines()

        elements = []
        atomic_num = []
        for line in lines[1:]:
            if line.split()[0] == 'ele':
                return elements, np.array(atomic_num, dtype=float)
            element = ''.join(re.findall(r'[a-zA-Z]', line.split()[0]))
            if element not in elements:
                elements.append(element)
                atomic_num.append(int(line.split()[1]))
        
        raise ValueError(
            f"Could not determine list of elements in file {self.wind_tablefile}.")

    def load_wind_yields(self) -> np.ndarray:
        """ Load tables of yields ejected as winds from massive stars.
        """

        wind_yld = np.zeros([len(self.elements), self.rot.size,
                             self.metal.size, self.mass.size])

        with open(self.wind_tablefile, 'r', encoding="utf-8") as file:
            lines = file.readlines()

        for index, line in enumerate(lines):
            if line.split()[0] == 'ele':
                model = line.split()[4]
                ind_metal = self._get_metal_index_from_model(model[3])
                ind_rot = self._get_rot_index_from_model(model[4:])

                data = np.genfromtxt(self.wind_tablefile,
                                     usecols=[1, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                                     skip_header=index+1,
                                     max_rows=142).T
                for i, atom_nr in enumerate(self.atomic_num):
                    mask = data[0] == atom_nr
                    wind_yld[i, ind_rot, ind_metal] = np.sum(
                        data[1:, mask], axis=1)

        return wind_yld

    def load_ccsn_yields(self) -> np.ndarray:
        """ Load tables of yields ejected by core-collapse supernovae.
        """

        wind_yld = self.load_wind_yields()
        total_yld = np.zeros([len(self.elements), self.rot.size,
                              self.metal.size, self.mass.size])

        with open(self.yield_tablefile, 'r', encoding="utf-8") as file:
            lines = file.readlines()

        for index, line in enumerate(lines):
            if line.split()[0] == 'ele':
                model = line.split()[4]
                ind_metal = self._get_metal_index_from_model(model[3])
                ind_rot = self._get_rot_index_from_model(model[4:])

                total_yld[:, ind_rot, ind_metal, :] = np.genfromtxt(
                    self.yield_tablefile,
                    usecols=[
                        4, 5, 6, 7, 8, 9, 10, 11, 12],
                    skip_header=index+1,
                    max_rows=53)

        ccsn_yld = total_yld - wind_yld
        ccsn_yld[ccsn_yld < 0.0] = 0.0
        ccsn_yld[:, :, :, (self.mass > self.ccsn_mmax)] = 0.0

        return ccsn_yld

    @staticmethod
    def _get_metal_index_from_model(model: str) -> int:
        """ Convenience function for converting table metal labels into table index.
        """
        if model == 'a':
            return 3
        elif model == 'b':
            return 2
        elif model == 'c':
            return 1
        elif model == 'd':
            return 0
        else:
            raise ValueError("Model does not exist.")

    @staticmethod
    def _get_rot_index_from_model(model: str) -> int:
        """ Convenience function for converting table rotation labels into table index.
        """
        if model == '000':
            return 0
        elif model == '150':
            return 1
        elif model == '300':
            return 2
        else:
            raise ValueError("Model does not exist.")
