"""
limongichieffi2018
==========

This submodule contains all the code required to load yields from
Limongi & Chieffi (2018).
"""

import io
import re
import tarfile
from typing import List

import astropy.units as u
import numpy as np
from astropy.units import Quantity

from ..population import StarPopulation
from .source import Source
from .yieldtables import YieldTables


class LimongiChieffi2018(YieldTables):
    """
    Include yields majority of elements (and isotopes) for massive stars
    [13, 15, 20, 25, 30, 40, 60, 80, 120] Msun, with initial iron abundance
    ([Fe/H]) [0, -1, -2, -3], and three different models for stellar
    rotation [0, 150, 300] km/s.

    Upon first initialization, the yield tables are automatically downloaded.
    If this fails, these files (.tgz) can be downloaded manually and placed/linked
    to directory LimongiChieffi2018 placed in the path of this file, e.g.,
    ln -s /path/to/downloaded/yields /path/to/arsenal_gear/element_yields/LimongiChieffi2018

    Yields available at https://orfeo.oa-roma.inaf.it/

    Reference: Limongi M & Chieffi A, 2018, ApJS, 237, 13L
    """

    # hard-coded parameters of the Limongi & Chieffi (2018) tables
    lc_url = "https://orfeo.oa-roma.inaf.it/"
    models = ["F", "I", "M", "R"]
    mass = np.array([13.0, 15.0, 20.0, 25.0, 30.0, 40.0, 60.0, 80.0, 120.0])
    metal = np.array([3.236e-5, 3.2363e-4, 3.236e-3, 1.345e-2])
    feh = np.array([-3, -2, -1, 0])
    rot = np.array([0, 150, 300])
    ccsn_mmax = 25

    def __init__(self, model: str = "R") -> None:
        """
        Args:
            model: choice of model to load, see Limongi & Chieffi (2018) for details

        Usage:
            >> lc2018 = arsenal_gear.element_yields.LimongiChieffi2018()
            >> mass = np.linspace(8, 120, 100)*u.M_sun
            >> metals = u.dimensionless_unscaled * 0.1 * np.ones(100)
            >> tform = u.Myr * np.zeros(100)
            >> rot = u.km / u.s * np.zeros(100)
            >> stars = arsenal_gear.population.StarPopulation(mass=mass, metals=metals, tform=tform, rot=rot)
            >> plt.plot(mass, yields.ccsn_yields('H', stars, interpolate='nearest'), '-')

        Attributes:
            lc_url           Yield source website
            model            Available models
            mass             Tabulated masses
            metal            Tabulated metallicities
            feh              Tabulated [Fe/H]
            rot              Tabulated stellar rotation velocities
            ccsn_max         Assumed minimum mass for direct collase to black hole
            filedir          Directory of this file (used for relative path)
            name             Name of yield tables
            elements         Elements available in table.
            atomic_num       Atomic numbers of available elements.
            wind             Source object for stellar winds (massive stars)
            ccsn             Source object for core-collapse SNe
        """

        if model not in self.models:
            raise ValueError(f"Model {model} does not exist.")

        super().__init__()
        self.filedir = self.filedir / "data/LimongiChieffi2018"
        self.name = "Limongi & Chieffi (2018)"

        if not self.filedir.is_dir():
            self.filedir.mkdir(parents=True, exist_ok=True)
            self.download_yields(table=f"tab_{model}")
        elif not (self.filedir / f"tab_{model}.tgz").is_file():
            self.download_yields(table=f"tab_{model}")

        self.filedir = self.filedir / f"tab_{model}.tgz"
        self.elements, self.atomic_num = self.get_element_list()

        # Stellar wind yields
        wind_yld = self.load_wind_yields()
        self.wind = Source(self.elements, [self.rot, self.metal, self.mass], wind_yld)

        # Core-collapse SNe yields
        ccsn_yld = self.load_ccsn_yields()
        self.ccsn = Source(self.elements, [self.rot, self.metal, self.mass], ccsn_yld)

    def ccsn_yields(
        self,
        elements: List[str],
        starPop: StarPopulation,
        interpolate: str = "nearest",
        extrapolate: bool = False,
    ) -> Quantity["mass"]:
        """Interpolate yields from core-collapse supernovae for specified elements.
        Stellar parameters can be provided as single value, array + single value, or arrays.

        Args:
            elements: list of elements, as specified by symbols (e.g., ['H'] for hydrogen)
            starSop: StarPopulation object
            interpolate: passed as method to scipy.interpolate.RegularGridInterpolator
            extrapolate: if False, then mass, metal, and rot are set to limits if outside bound
        Returns:
            List of yields matching provided element list

        """
        args = [
            starPop["rot"].to(u.km / u.s).value,
            starPop["metals"].value,
            starPop["mass"].to(u.M_sun).value,
        ]
        return (
            self.ccsn.get_yld(
                elements, args, interpolate=interpolate, extrapolate=extrapolate
            )
            * u.M_sun
        )

    def wind_yields(
        self,
        elements: List[str],
        starPop: StarPopulation,
        interpolate: str = "nearest",
        extrapolate: bool = False,
    ) -> Quantity["mass"]:
        """Interpolate yields from massive stars ejected as winds for specified elements.
        Stellar parameters can be provided as single value, array + single value, or arrays.

        Args:
            elements: list of elements, as specified by symbols (e.g., ['H'] for hydrogen)
            starSop: StarPopulation object
            interpolate: passed as method to scipy.interpolate.RegularGridInterpolator
            extrapolate: if False, then mass, metal, and rot are set to table if outside bound.
        Returns:
            List of yields matching provided element list

        """
        args = [
            starPop["rot"].to(u.km / u.s).value,
            starPop["metals"].value,
            starPop["mass"].to(u.M_sun).value,
        ]
        return (
            self.wind.get_yld(
                elements, args, interpolate=interpolate, extrapolate=extrapolate
            )
            * u.M_sun
        )

    def get_element_list(self) -> None:
        """Read element symbols and atomic numbers from tables."""

        elements = []
        atomic_num = []

        with tarfile.open(self.filedir, "r:*") as tar:
            windfile = tar.extractfile("tabwind.dec")
            lines = io.TextIOWrapper(windfile, encoding="utf-8").readlines()

        for line in lines[1:]:
            if line.split()[0] == "ele":
                return elements, np.array(atomic_num, dtype=float)
            element = "".join(re.findall(r"[a-zA-Z]", line.split()[0]))
            if element not in elements:
                elements.append(element)
                atomic_num.append(int(line.split()[1]))

        raise ValueError("Could not determine list of elements from tabwind.dec")

    def load_wind_yields(self) -> np.ndarray:
        """Load tables of yields ejected as winds from massive stars."""

        wind_yld = np.zeros(
            [len(self.elements), self.rot.size, self.metal.size, self.mass.size]
        )

        with tarfile.open(self.filedir, "r:*") as tar:
            f = tar.extractfile("tabwind.dec")
            lines = io.TextIOWrapper(f, encoding="utf-8").readlines()

        for index, line in enumerate(lines):
            if line.split()[0] == "ele":
                model = line.split()[4]
                ind_metal = self._get_metal_index_from_model(model[3])
                ind_rot = self._get_rot_index_from_model(model[4:])

                with tarfile.open(self.filedir, "r:*") as tar:
                    f = tar.extractfile("tabwind.dec")
                    data = np.genfromtxt(
                        f,
                        usecols=[1, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                        skip_header=index + 1,
                        max_rows=142,
                    ).T

                for i, atom_nr in enumerate(self.atomic_num):
                    mask = data[0] == atom_nr
                    wind_yld[i, ind_rot, ind_metal] = np.sum(data[1:, mask], axis=1)

        return wind_yld

    def load_ccsn_yields(self) -> np.ndarray:
        """Load tables of yields ejected by core-collapse supernovae."""

        wind_yld = self.load_wind_yields()
        total_yld = np.zeros(
            [len(self.elements), self.rot.size, self.metal.size, self.mass.size]
        )

        with tarfile.open(self.filedir, "r:*") as tar:
            f = tar.extractfile("tab_yieldstot_ele_exp.dec")
            lines = io.TextIOWrapper(f, encoding="utf-8").readlines()

        for index, line in enumerate(lines):
            if line.split()[0] == "ele":
                model = line.split()[4]
                ind_metal = self._get_metal_index_from_model(model[3])
                ind_rot = self._get_rot_index_from_model(model[4:])

                with tarfile.open(self.filedir, "r:*") as tar:
                    f = tar.extractfile("tab_yieldstot_ele_exp.dec")
                    total_yld[:, ind_rot, ind_metal, :] = np.genfromtxt(
                        f,
                        usecols=[4, 5, 6, 7, 8, 9, 10, 11, 12],
                        skip_header=index + 1,
                        max_rows=53,
                    )

        ccsn_yld = total_yld - wind_yld
        ccsn_yld[ccsn_yld < 0.0] = 0.0
        ccsn_yld[:, :, :, (self.mass > self.ccsn_mmax)] = 0.0

        return ccsn_yld

    def download_yields(self, table: str) -> None:
        """Downloader for tabulated yield files."""
        from ..utils.scraper import downloader

        downloader(
            self.filedir / f"{table}.tgz",
            f"{self.lc_url}/2018-modelli/yields/{table}.tgz",
            message=f"Yield file {table}.tgz not found.",
        )

        if not (self.filedir / "readme.txt").is_file():
            downloader(
                self.filedir / "readme.txt",
                f"{self.lc_url}/2018-modelli/yields/legenda",
                message="See downloaded readme.txt for info about yields.",
            )

    @staticmethod
    def _get_metal_index_from_model(model: str) -> int:
        """Convenience function for converting table metal labels into table index."""
        met_ind_mods = {"a": 3, "b": 2, "c": 1, "d": 0}
        if model in met_ind_mods:
            return met_ind_mods[model]
        else:
            raise ValueError("Model does not exist.")

    @staticmethod
    def _get_rot_index_from_model(model: str) -> int:
        """Convenience function for converting table rotation labels into table index."""
        rot_ind_mods = {"000": 0, "150": 1, "300": 2}
        if model in rot_ind_mods:
            return rot_ind_mods[model]
        else:
            raise ValueError("Model does not exist.")
