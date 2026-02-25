"""
doherty
==========

This submodule contains all the code required to load yields from
Doherty et al. (2014a,b).
"""

import re
import zipfile
from typing import List

import astropy.units as u
import numpy as np
from astropy.units import Quantity

from ..population import StarPopulation
from .source import Source
from .yieldtables import YieldTables


class Doherty2014a(YieldTables):
    """
    Yields available in suplementary material of Doherty et al. (2014a)

    Reference: Doherty et al., 2014a, MNRAS, 437, 195
    """

    models = ["VW93", "VW-M"]
    d14_url = (
        "https://academic.oup.com/mnras/article/441/1/582/980303#supplementary-data"
    )
    file = "998013_Supplementary_Data.zip"
    mass = np.array([6.5, 7.0, 7.5, 8.0, 8.5, 9.0])
    metal = np.array([0.004, 0.008, 0.02])

    def __init__(self, model: str = "VW93"):
        """
        Args:
            model: Model used for mass loss rates.
        Usage:
            >> d14a = arsenal_gear.element_yields.Doherty2014a()
            >> mass = np.linspace(8, 120, 100)*u.M_sun
            >> metals = u.dimensionless_unscaled * 0.1 * np.ones(100)
            >> tform = u.Myr * np.zeros(100)
            >> stars = arsenal_gear.population.StarPopulation(mass=mass, metals=metals, tform=tform)
            >> plt.plot(mass, d14a.ccsn_yields('H', stars, interpolate='nearest'), '-')

        Attributes:
            d14_url          Yield source website
            mass             Tabulated masses
            metal            Tabulated metallicities
            filedir          Directory of yield tables
            name             Name of yields tables
            elements         Elements available in table
            atomic_num       Atomic numbers of available elements
            sagb             Source object for stellar winds (super asymptotic giant branch)
        """
        if model not in self.models:
            raise ValueError(f"Model {model} does not exist.")
        self.model = model

        super().__init__()
        self.filedir = self.filedir / "data" / "Doherty"
        self.name = "Doherty et al. (2014a)"
        self.filename = "TABLE1-VW93ML.txt"
        if self.model == "VW-M":
            self.filename = "TABLE2-VW-MML.txt"

        if not self.filedir.is_dir():
            self.filedir.mkdir(parents=True, exist_ok=True)
            self.download_yields()
        else:
            if not (self.filedir / self.file).is_file():
                self.download_yields()
        self.filedir = self.filedir / self.file

        self.elements, self.atomic_num = self.get_element_list()

        # Super asymptotic giant branch
        sagb_yld = self.load_sagb_yields()
        self.sagb = Source(self.elements, [self.metal, self.mass], sagb_yld)

    def sagb_yields(
        self,
        elements: List[str],
        starPop: StarPopulation,
        interpolate: str = "nearest",
        extrapolate: bool = False,
    ) -> Quantity["mass"]:
        """Interpolate yields from stars on the super asymptotic giant branch for specified elements.
        Stellar parameters can be provided as single value, array + single value, or arrays.

        Args:
            elements: list of elements, as specified by symbols (e.g., ['H'] for hydrogen)
            starPop: stellar population object
            interpolate: passed as method to scipy.interpolate.RegularGridInterpolator
            extrapolate: if False, then mass, metal, and rot are set to table if outside bound
        Returns:
            Dictionary of yields matching provided element list

        """
        elements = np.atleast_1d(elements)
        args = [
            starPop["metals"].value,
            starPop["mass"].to(u.M_sun).value,
        ]
        yld_array = self.sagb.get_yld(
            elements, args, interpolate=interpolate, extrapolate=extrapolate
        )

        if len(elements) == 1:
            return {elements[0]: yld_array * u.M_sun}
        return {element: yld_array[i] * u.M_sun for i, element in enumerate(elements)}

    ## Internal functions for loading data
    def get_element_list(self) -> None:
        """Read element symbols and atomic numbers from tables."""

        elements = []
        atomic_num = None

        with zipfile.ZipFile(self.filedir) as zipf:
            with zipf.open(self.filename) as sagbfile:
                lines = sagbfile.read().decode("utf-8").splitlines()

        for line in lines[2:]:
            if line.split()[0] == "g":
                return elements, atomic_num
            if line.split()[0] == "p":
                element = "H"
            else:
                element = "".join(re.findall(r"[a-zA-Z]", line.split()[0]))
                element = element[:1].upper() + element[1:]
            if element not in elements:
                elements.append(element)

        raise ValueError("Could not determine list of elements from {self.filename}")

    def load_sagb_yields(self) -> np.ndarray:
        """Load tables of yields ejected as winds from super asymptotic giant branch stars."""

        sagb_yld = np.zeros([len(self.elements), self.metal.size, self.mass.size])

        with zipfile.ZipFile(self.filedir) as zipf:
            with zipf.open(self.filename) as sagbfile:
                lines = sagbfile.read().decode("utf-8").splitlines()

                ind_mass = None
                ind_metal = None
                for line in lines:
                    if "alpha" in line:
                        continue
                    if self.model in line:
                        mass = float(line.split("M")[0])
                        metal = float(line.split("Z=")[-1].split(self.model)[0])
                        ind_mass = np.argmin(np.abs(self.mass - mass))
                        ind_metal = np.argmin(np.abs(self.metal - metal))
                        continue
                    if any(s in line for s in ("Species", "g", "d")):
                        continue
                    if line.split()[0] == "p":
                        ind_elem = 0
                    else:
                        current_elem = "".join(re.findall(r"[a-zA-Z]", line.split()[0]))
                        current_elem = current_elem[:1].upper() + current_elem[1:]
                        ind_elem = next(
                            i
                            for i, elem in enumerate(self.elements)
                            if current_elem == elem
                        )

                    sagb_yld[ind_elem, ind_metal, ind_mass] += float(line.split()[2])

        return sagb_yld

    def download_yields(self) -> None:
        """Downloader for tabulated yield files."""
        # from ..utils.scraper import downloader

        print("Doherty (2014a) tables not found, starting download.")
        raise RuntimeError(
            f"Automatic download for these yields is not available. "
            f"Please download yields from suplementary material in {self.name} and place"
            f"in {self.filedir} (follow link {self.d14_url} for download)"
        )
        # downloader(self.filedir / self.file, self.d14_url, message=None)


class Doherty2014b(YieldTables):
    """
    Yields available in suplementary material of Doherty et al. (2014b)

    Reference: Doherty et al., 2014b, MNRAS, 441, 582
    """

    models = ["VW93", "B95"]
    d14_url = (
        "https://academic.oup.com/mnras/article/441/1/582/980303#supplementary-data"
    )
    file = "stu571_Supplementary_Data.zip"
    mass = np.array([6.5, 7.0, 7.5])
    metal = np.array([0.0001, 0.001])

    def __init__(self, model: str = "VW93"):
        """
        Args:
            model: Model used for mass loss rates.
        Usage:
            >> d14b = arsenal_gear.element_yields.Doherty2014b()
            >> mass = np.linspace(8, 120, 100)*u.M_sun
            >> metals = u.dimensionless_unscaled * 0.1 * np.ones(100)
            >> tform = u.Myr * np.zeros(100)
            >> stars = arsenal_gear.population.StarPopulation(mass=mass, metals=metals, tform=tform)
            >> plt.plot(mass, d14b.ccsn_yields('H', stars, interpolate='nearest'), '-')

        Attributes:
            d14_url          Yield source website
            mass             Tabulated masses
            metal            Tabulated metallicities
            filedir          Directory of yield tables
            name             Name of yields tables
            elements         Elements available in table
            atomic_num       Atomic numbers of available elements
            sagb             Source object for stellar winds (super asymptotic giant branch)
        """
        if model not in self.models:
            raise ValueError(f"Model {model} does not exist.")
        self.model = model

        super().__init__()
        self.filedir = self.filedir / "data" / "Doherty"
        self.name = "Doherty et al. (2014b)"
        self.filename = "P3Doh14b-table1.txt"
        if self.model == "B95":
            self.filename = "P3Doh14b-table2.txt"

        if not self.filedir.is_dir():
            self.filedir.mkdir(parents=True, exist_ok=True)
            self.download_yields()
        else:
            if not (self.filedir / self.file).is_file():
                self.download_yields()
        self.filedir = self.filedir / self.file

        self.elements, self.atomic_num = self.get_element_list()

        # Super asymptotic giant branch
        sagb_yld = self.load_sagb_yields()
        self.sagb = Source(self.elements, [self.metal, self.mass], sagb_yld)

    def sagb_yields(
        self,
        elements: List[str],
        starPop: StarPopulation,
        interpolate: str = "nearest",
        extrapolate: bool = False,
    ) -> Quantity["mass"]:
        """Interpolate yields from stars on the super asymptotic giant branch for specified elements.
        Stellar parameters can be provided as single value, array + single value, or arrays.

        Args:
            elements: list of elements, as specified by symbols (e.g., ['H'] for hydrogen)
            starPop: stellar population object
            interpolate: passed as method to scipy.interpolate.RegularGridInterpolator
            extrapolate: if False, then mass, metal, and rot are set to table if outside bound
        Returns:
            Dictionary of yields matching provided element list

        """
        elements = np.atleast_1d(elements)
        args = [
            starPop["metals"].value,
            starPop["mass"].to(u.M_sun).value,
        ]
        yld_array = self.sagb.get_yld(
            elements, args, interpolate=interpolate, extrapolate=extrapolate
        )

        if len(elements) == 1:
            return {elements[0]: yld_array * u.M_sun}
        return {element: yld_array[i] * u.M_sun for i, element in enumerate(elements)}

    ## Internal functions for loading data
    def get_element_list(self) -> None:
        """Read element symbols and atomic numbers from tables."""
        return Doherty2014a.get_element_list(self)

    def load_sagb_yields(self) -> np.ndarray:
        """Load tables of yields ejected as winds from super asymptotic giant branch stars."""
        return Doherty2014a.load_sagb_yields(self)

    def download_yields(self) -> None:
        """Downloader for tabulated yield files."""
        # from ..utils.scraper import downloader

        print("Doherty (2014b) tables not found, starting download.")
        raise RuntimeError(
            f"Automatic download for these yields is not available. "
            f"Please download yields from suplementary material in {self.name} and place"
            f"in {self.filedir} (follow link {self.d14_url} for download)"
        )
        # downloader(self.filedir / self.file, self.d14_url, message=None)
