"""
nugrid
==========

This submodule contains all the code required to load yields from
NuGrid.
"""

import os
import re
from typing import List

import astropy.units as u
import numpy as np
from astropy.units import Quantity

from ..population import StarPopulation
from .source import Source
from .yieldtables import YieldTables


class NuGrid(YieldTables):
    """
    NuGrid is a collection of yields, creating a regular grid for masses
    [1.0, 1.65, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 12.0, 15.0, 20.0, 25.0] Msun
    and metallicities [0.0001, 0.001, 0.006, 0.01, 0.02] by combining yields
    from Ritter et al. (2018; see also, Pignatari et al., 2016) and
    Battino et al. (2019, 2021).
    See classes Ritter2018 and Battino20192021 for more details.

    Yields available at https://nugrid.github.io/

    References: Pignatari et al., 2016, ApJS, 225, 24; Ritter et al., 2018, MNRAS, 480, 538;
                Battino et al., 2019, MNRAS, 489, 1082; Battino et al., 2021, Universe, 7, 25;
    """

    # hard-coded parameters of the NuGrid tables
    ng_url = "https://nugrid.github.io/content/data/"
    models = ["delay", "rapid"]
    lo_mass = np.array([1.0, 1.65, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
    hi_mass = np.array([12.0, 15.0, 20.0, 25.0])
    metal = np.array([0.0001, 0.001, 0.006, 0.01, 0.02])

    def __init__(self, model: str = "delay") -> None:
        """
        Args:
            model: choise of model to load, see Fryer et al. (2012) for details.

        Usage:
            >> nugrid = arsenal_gear.element_yields.NuGrid()
            >> mass = np.linspace(8, 120, 1000)*u.M_sun
            >> metals = u.dimensionless_unscaled * 0.1 * np.ones(100)
            >> tform = u.Myr * np.zeros(100)
            >> stars = arsenal_gear.population.StarPopulation(mass=mass, metals=metals, tform=tform)
            >> plt.plot(mass, nugrid.ccsn_yields('H', stars,
                                                 interpolate='nearest'),
                        '-', color=colors[-1])

        Attributes:
            models           Available models.
            mass             Tabulated masses.
            metal            Tabulated metallicities.
            filedir          Directory of this file (used for relative path).
            elements         Elements available in table.
            atomic_num       Atomic numbers of available elements.
            agb              Source object for stellar winds (asymptotic giant branch)
            wind             Source object for stellar winds (massive stars)
            ccsn             Source object for core-collapse SNe
        """
        if model not in self.models:
            raise ValueError(f"Model {model} does not exist.")

        super().__init__()
        self.filedir += "/NuGrid"
        self.name = "NuGrid"

        r18 = Ritter2018(model=model)
        bat = Battino20192021()

        self.elements = r18.elements
        self.atomic_num = r18.elements

        # Asymptotic giant branch star yields
        agb_yld = r18.load_agb_yields()
        agb_yld_bat = bat.load_agb_yields()
        bat_idx = np.array(
            [np.where(np.asarray(bat.elements) == el)[0][0] for el in r18.elements]
        )

        r18_replace = np.array([(1, 2), (3, 2), (4, 2), (1, 3), (3, 3), (4, 3)])
        bat_replace = np.array([(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1)])
        ri, rj = r18_replace.T
        bi, bj = bat_replace.T
        agb_yld[:, ri, rj] = agb_yld_bat[bat_idx[:, None], bi, bj]
        self.agb = Source(self.elements, [self.metal, self.lo_mass], agb_yld)

        # Stellar wind yields
        self.wind = r18.wind

        # Core-collapse SNe yields
        self.ccsn = r18.ccsn

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
            elements: list of elements, as specified by symbols (e.g., ['H'] for hydrogen).
            starPop: StarPopulation object
            interpolate: passed as method to scipy.interpolate.RegularGridInterpolator
            extrapolate: if False, then mass, metal, and rot are set to limits if outside bound.
        Returns:
            List of yields matching provided element list

        """
        args = [
            starPop["metals"].value,
            starPop["mass"].to(u.M_sun).value,
        ]

        yields = np.zeros_like(starPop["mass"].value)
        yields = self.ccsn.get_yld(
            elements, args, interpolate=interpolate, extrapolate=extrapolate
        )
        return yields * u.M_sun

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
            elements: list of elements, as specified by symbols (e.g., ['H'] for hydrogen).
            starPop: StarPopulation object
            interpolate: passed as method to scipy.interpolate.RegularGridInterpolator
            extrapolate: if False, then mass, metal, and rot are set to table if outside bound.
        Returns:
            List of yields matching provided element list

        """
        args = [
            starPop["metals"].value,
            starPop["mass"].to(u.M_sun).value,
        ]

        yields = np.zeros_like(starPop["mass"].value)
        yields = self.wind.get_yld(
            elements, args, interpolate=interpolate, extrapolate=extrapolate
        )
        return yields * u.M_sun

    def agb_yields(
        self,
        elements: List[str],
        starPop: StarPopulation,
        interpolate: str = "nearest",
        extrapolate: bool = False,
    ) -> Quantity["mass"]:
        """Interpolate yields from AGB stars ejected as winds for specified elements.
        Stellar parameters can be provided as single value, array + single value, or arrays.

        Args:
            elements: list of elements, as specified by symbols (e.g., ['H'] for hydrogen).
            starPop: StarPopulation object
            interpolate: passed as method to scipy.interpolate.RegularGridInterpolator
            extrapolate: if False, then mass, metal, and rot are set to table if outside bound.
        Returns:
            List of yields matching provided element list

        """
        args = [
            starPop["metals"].value,
            starPop["mass"].to(u.M_sun).value,
        ]

        yields = np.zeros_like(starPop["mass"].value)
        yields = self.agb.get_yld(
            elements, args, interpolate=interpolate, extrapolate=extrapolate
        )
        return yields * u.M_sun


class Pignatari2016(YieldTables):
    """
    Include yields for many elements for stars with masses
    [1.65, 2.0, 3.0, 5.0, 15, 20, 25] Msun, with initial metallicities
    [0.01, 0.02]. All explosion models assumes delayed trigger (see Fryer et al., 2012).

    Upon first initialization, the yield tables are automatically downloaded.
    If this fails, these files (.txt) can be downloaded manually and placed/linked
    to directory NuGrid placed in the path of this file, e.g.,
    ln -s /path/to/downloaded/yields /path/to/arsenal_gear/element_yields/NuGrid

    Yields available at https://nugrid.github.io/

    References: Pignatari et al., 2016, ApJS, 225, 24
    """

    # hard-coded parameters of the Pignatari (2016) tables
    p16_url = "https://download.nugridstars.org/set1/Yield_tables/"
    files = [
        "element_table_set1.1_yields_exp.txt",
        "element_table_set1.1_yields_winds.txt",
        "element_table_set1.2_yields_exp.txt",
        "element_table_set1.2_yields_winds.txt",
    ]
    lo_mass = np.array([1.65, 2.0, 3.0, 5.0])
    hi_mass = np.array([15.0, 20.0, 25.0])
    metal = np.array([0.01, 0.02])

    def __init__(self) -> None:

        super().__init__()
        self.filedir += "/NuGrid"
        self.name = "Pignatari et al. (2016)"

        if not os.path.isdir(self.filedir):
            os.mkdir(self.filedir)
            self.download_yields()
        else:
            if not all(
                [os.path.isfile(self.filedir + os.sep + file) for file in self.files]
            ):
                self.download_yields()

        self.elements, self.atomic_num = self.get_element_list()

        # Asymptotic giant branch
        agb_yld = self.load_agb_yields()
        self.agb = Source(self.elements, [self.metal, self.lo_mass], agb_yld)

        # Winds (massive stars)
        wind_yld = self.load_wind_yields()
        self.wind = Source(self.elements, [self.metal, self.hi_mass], wind_yld)

        # Core-collapse SNe
        ccsn_yld = self.load_ccsn_yields()
        self.ccsn = Source(self.elements, [self.metal, self.hi_mass], ccsn_yld)

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
            elements: list of elements, as specified by symbols (e.g., ['H'] for hydrogen).
            starPop: stellar population object
            interpolate: passed as method to scipy.interpolate.RegularGridInterpolator
            extrapolate: if False, then mass, metal, and rot are set to limits if outside bound.
        Returns:
            List of yields matching provided element list

        """
        args = [
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
            elements: list of elements, as specified by symbols (e.g., ['H'] for hydrogen).
            starPop: stellar population object
            interpolate: passed as method to scipy.interpolate.RegularGridInterpolator
            extrapolate: if False, then mass, metal, and rot are set to table if outside bound.
        Returns:
            List of yields matching provided element list

        """
        args = [
            starPop["metals"].value,
            starPop["mass"].to(u.M_sun).value,
        ]
        return (
            self.wind.get_yld(
                elements, args, interpolate=interpolate, extrapolate=extrapolate
            )
            * u.M_sun
        )

    def agb_yields(
        self,
        elements: List[str],
        starPop: StarPopulation,
        interpolate: str = "nearest",
        extrapolate: bool = False,
    ) -> Quantity["mass"]:
        """Interpolate yields from stars on the asymptotic giant branch for specified elements.
        Stellar parameters can be provided as single value, array + single value, or arrays.

        Args:
            elements: list of elements, as specified by symbols (e.g., ['H'] for hydrogen).
            starPop: stellar population object
            interpolate: passed as method to scipy.interpolate.RegularGridInterpolator
            extrapolate: if False, then mass, metal, and rot are set to table if outside bound.
        Returns:
            List of yields matching provided element list

        """
        args = [
            starPop["metals"].value,
            starPop["mass"].to(u.M_sun).value,
        ]
        return (
            self.agb.get_yld(
                elements, args, interpolate=interpolate, extrapolate=extrapolate
            )
            * u.M_sun
        )

    def get_element_list(self) -> None:
        """Read element symbols from tables, atomic numbers not available."""

        elements = []
        atomic_numbers = []

        with open(self.filedir + os.sep + self.files[0], encoding="utf-8") as f:
            lines = f.readlines()

        for line in lines[3:]:
            element = "".join(re.findall(r"[a-zA-Z]", line.split("&")[1]))
            if element not in elements:
                elements.append(element)

        return elements, atomic_numbers

    def load_agb_yields(self) -> np.ndarray:
        """Load tables of yields ejected during agb."""

        agb_yld = np.zeros([len(self.elements), self.metal.size, self.lo_mass.size])

        for file in self.files:
            if not "winds" in file:
                continue

            if "set1.1" in file:
                ind_met = 0
            elif "set1.2" in file:
                ind_met = 1
            else:
                raise ValueError(
                    "Could not find file while loading Pignatari et al. (2016) yields"
                )

            data = np.genfromtxt(
                self.filedir + os.sep + file,
                usecols=[2, 3, 4, 5],
                skip_header=3,
                delimiter="&",
            )

            agb_yld[:, ind_met, :] = data

        return agb_yld

    def load_wind_yields(self) -> np.ndarray:
        """Load tables of yields ejected as winds from massive stars."""

        wind_yld = np.zeros([len(self.elements), self.metal.size, self.hi_mass.size])

        for file in self.files:
            if not "winds" in file:
                continue

            if "set1.1" in file:
                ind_met = 0
            elif "set1.2" in file:
                ind_met = 1
            else:
                raise ValueError(
                    "Could not find file while loading Pignatari et al. (2016) yields"
                )

            data = np.genfromtxt(
                self.filedir + os.sep + file,
                usecols=[6, 7, 8],
                skip_header=3,
                delimiter="&",
            )

            wind_yld[:, ind_met, :] = data

        return wind_yld

    def load_ccsn_yields(self) -> np.ndarray:
        """Load tables of yields ejected by core-collapse supernovae."""

        ccsn_yld = np.zeros([len(self.elements), self.metal.size, self.hi_mass.size])

        for file in self.files:
            if not "exp" in file:
                continue

            if "set1.1" in file:
                ind_met = 0
            elif "set1.2" in file:
                ind_met = 1
            else:
                raise ValueError(
                    "Could not find file while loading Pignatari et al. (2016) yields"
                )

            data = np.genfromtxt(
                self.filedir + os.sep + file,
                usecols=[2, 4, 6],
                skip_header=3,
                delimiter="&",
            )

            ccsn_yld[:, ind_met, :] = data

        return ccsn_yld

    def download_yields(self) -> None:
        """Downloader for tabulated yield files."""
        from ..utils.scraper import downloader

        print("NuGrid (Pignatari et al., 2016) tables not found, starting download.")
        for file in self.files:
            downloader(
                self.filedir + os.sep + file, self.p16_url + os.sep + file, message=None
            )


class Ritter2018(YieldTables):
    """
    Include yields for many elements for stars with masses
    [1.0, 1.65, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 12, 15, 20, 25] Msun, with initial metallicities
    [0.0001, 0.001, 0.006, 0.01, 0.02].
    Includes option for both rapid and delayed explosion trigger (see Fryer et al., 2012).

    Upon first initialization, the yield tables are automatically downloaded.
    If this fails, these files (.txt) can be downloaded manually and placed/linked
    to directory NuGrid placed in the path of this file, e.g.,
    ln -s /path/to/downloaded/yields /path/to/arsenal_gear/element_yields/NuGrid

    Yields available at https://nugrid.github.io/

    References:  Ritter et al., 2018, MNRAS, 480, 538
    """

    # hard-coded parameters of the Ritter (2018) tables
    r18_url = "https://download.nugridstars.org/set1ext/Yield_tables/"
    files = [
        "element_yield_table_MESAonly_fryer12_delay_total.txt",
        "element_yield_table_MESAonly_fryer12_delay_winds.txt",
        "element_yield_table_MESAonly_fryer12_rapid_total.txt",
        "element_yield_table_MESAonly_fryer12_rapid_winds.txt",
    ]
    models = ["delay", "rapid"]
    lo_mass = np.array([1.0, 1.65, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
    hi_mass = np.array([12.0, 15.0, 20.0, 25.0])
    metal = np.array([0.0001, 0.001, 0.006, 0.01, 0.02])

    def __init__(self, model: str = "delay") -> None:

        if model not in self.models:
            raise ValueError(f"Model {model} does not exist.")

        self.files = [file for file in self.files if model in file]

        super().__init__()
        self.filedir += "/NuGrid"
        self.name = "Ritter et al. (2018)"

        if not os.path.isdir(self.filedir):
            os.mkdir(self.filedir)
            self.download_yields()
        else:
            if not all(
                [os.path.isfile(self.filedir + os.sep + file) for file in self.files]
            ):
                self.download_yields()

        self.elements, self.atomic_num = self.get_element_list()

        # Asymptotic giant branch
        agb_yld = self.load_agb_yields()
        self.agb = Source(self.elements, [self.metal, self.lo_mass], agb_yld)

        # Winds (massive stars)
        wind_yld = self.load_wind_yields()
        self.wind = Source(self.elements, [self.metal, self.hi_mass], wind_yld)

        # Core-collapse SNe
        ccsn_yld = self.load_ccsn_yields()
        self.ccsn = Source(self.elements, [self.metal, self.hi_mass], ccsn_yld)

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
            elements: list of elements, as specified by symbols (e.g., ['H'] for hydrogen).
            starPop: stellar population object
            interpolate: passed as method to scipy.interpolate.RegularGridInterpolator
            extrapolate: if False, then mass, metal, and rot are set to limits if outside bound.
        Returns:
            List of yields matching provided element list

        """
        args = [
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
            elements: list of elements, as specified by symbols (e.g., ['H'] for hydrogen).
            starPop: stellar population object
            interpolate: passed as method to scipy.interpolate.RegularGridInterpolator
            extrapolate: if False, then mass, metal, and rot are set to table if outside bound.
        Returns:
            List of yields matching provided element list

        """
        args = [
            starPop["metals"].value,
            starPop["mass"].to(u.M_sun).value,
        ]
        return (
            self.wind.get_yld(
                elements, args, interpolate=interpolate, extrapolate=extrapolate
            )
            * u.M_sun
        )

    def agb_yields(
        self,
        elements: List[str],
        starPop: StarPopulation,
        interpolate: str = "nearest",
        extrapolate: bool = False,
    ) -> Quantity["mass"]:
        """Interpolate yields from stars on the asymptotic giant branch for specified elements.
        Stellar parameters can be provided as single value, array + single value, or arrays.

        Args:
            elements: list of elements, as specified by symbols (e.g., ['H'] for hydrogen).
            starPop: stellar population object
            interpolate: passed as method to scipy.interpolate.RegularGridInterpolator
            extrapolate: if False, then mass, metal, and rot are set to table if outside bound.
        Returns:
            List of yields matching provided element list

        """
        args = [
            starPop["metals"].value,
            starPop["mass"].to(u.M_sun).value,
        ]
        return (
            self.agb.get_yld(
                elements, args, interpolate=interpolate, extrapolate=extrapolate
            )
            * u.M_sun
        )

    def get_element_list(self) -> None:
        """Read element symbols and atomic numbers from tables."""

        elements = []
        atomic_numbers = []

        with open(self.filedir + os.sep + self.files[0], encoding="utf-8") as f:
            lines = f.readlines()

        for line in lines[10:]:
            if "Table" in line:
                return elements, atomic_numbers
            element = "".join(re.findall(r"[a-zA-Z]", line.split("&")[1]))
            if element not in elements:
                elements.append(element)
                atomic_numbers.append(int(line.split("&")[-1]))

    def load_agb_yields(self) -> np.ndarray:
        """Load tables of yields ejected during agb."""

        agb_yld = np.zeros([len(self.elements), self.metal.size, self.lo_mass.size])

        for file in self.files:
            if not "winds" in file:
                continue

            with open(self.filedir + os.sep + file, encoding="utf-8") as f:
                lines = f.readlines()

            for iline, line in enumerate(lines):
                if "H Table" in line:
                    m = float(line.split()[-1].split(",")[0].split("=")[-1])
                    Z = float(line.split()[-1].split(",")[-1][:-1].split("=")[-1])
                    if m > max(self.lo_mass):
                        continue
                    ind_mass = np.argwhere(self.lo_mass == m)[0][0]
                    ind_met = np.argwhere(self.metal == Z)[0][0]

                    data = np.genfromtxt(
                        self.filedir + os.sep + file,
                        usecols=[2],
                        skip_header=iline + 4,
                        max_rows=len(self.elements),
                        delimiter="&",
                    )
                    agb_yld[:, ind_met, ind_mass] = data

        return agb_yld

    def load_wind_yields(self) -> np.ndarray:
        """Load tables of yields ejected as winds from massive stars."""

        wind_yld = np.zeros([len(self.elements), self.metal.size, self.hi_mass.size])

        for file in self.files:
            if not "winds" in file:
                continue

            with open(self.filedir + os.sep + file, encoding="utf-8") as f:
                lines = f.readlines()

            for iline, line in enumerate(lines):
                if "H Table" in line:
                    m = float(line.split()[-1].split(",")[0].split("=")[-1])
                    Z = float(line.split()[-1].split(",")[-1][:-1].split("=")[-1])
                    if m < min(self.hi_mass):
                        continue
                    ind_mass = np.argwhere(self.hi_mass == m)[0][0]
                    ind_met = np.argwhere(self.metal == Z)[0][0]

                    data = np.genfromtxt(
                        self.filedir + os.sep + file,
                        usecols=[2],
                        skip_header=iline + 4,
                        max_rows=len(self.elements),
                        delimiter="&",
                    )
                    wind_yld[:, ind_met, ind_mass] = data

        return wind_yld

    def load_ccsn_yields(self) -> np.ndarray:
        """Load tables of yields ejected by core-collapse supernovae."""

        total_yld = np.zeros([len(self.elements), self.metal.size, self.hi_mass.size])

        for file in self.files:
            if not "total" in file:
                continue

            with open(self.filedir + os.sep + file, encoding="utf-8") as f:
                lines = f.readlines()

            for iline, line in enumerate(lines):
                if "H Table" in line:
                    m = float(line.split()[-1].split(",")[0].split("=")[-1])
                    Z = float(line.split()[-1].split(",")[-1][:-1].split("=")[-1])
                    if m < min(self.hi_mass):
                        continue
                    ind_mass = np.argwhere(self.hi_mass == m)[0][0]
                    ind_met = np.argwhere(self.metal == Z)[0][0]

                    data = np.genfromtxt(
                        self.filedir + os.sep + file,
                        usecols=[2],
                        skip_header=iline + 4,
                        max_rows=len(self.elements),
                        delimiter="&",
                    )
                    total_yld[:, ind_met, ind_mass] = data
        wind_yld = self.load_wind_yields()

        return total_yld - wind_yld

    def download_yields(self) -> None:
        """Downloader for tabulated yield files."""
        from ..utils.scraper import downloader

        print("NuGrid (Ritter et al., 2018) tables not found, starting download.")
        for file in self.files:
            downloader(
                self.filedir + os.sep + file, self.r18_url + os.sep + file, message=None
            )


class Battino20192021(YieldTables):
    """
    Include yields for many elements for stars with masses [2.0, 3.0] Msun,
    with initial metallicities [0.001, 0.01, 0.02, 0.03].

    Upon first initialization, the yield tables are automatically downloaded.
    If this fails, these files (.txt) can be downloaded manually and placed/linked
    to directory NuGrid placed in the path of this file, e.g.,
    ln -s /path/to/downloaded/yields /path/to/arsenal_gear/element_yields/NuGrid

    Yields available at https://nugrid.github.io/

    References: Battino et al., 2019, MNRAS, 489, 1082; Battino et al., 2021, Universe, 7, 25
    """

    # hard-coded parameters of the Battino (2019,2021) tables
    b19_url = (
        "https://download.nugridstars.org/set1upd/yields_finalabu_tables/Battino2019"
    )
    b21_url = (
        "https://download.nugridstars.org/set1upd/yields_finalabu_tables/Battino2021"
    )
    files = ["Battino2019_Yields_table.txt", "Battino2021_Yields_table.txt"]
    mass = np.array([2.0, 3.0])
    metal = np.array([0.001, 0.01, 0.02, 0.03])

    def __init__(self) -> None:

        super().__init__()
        self.filedir += "/NuGrid"
        self.name = "Battino et al. (2019, 2021)"

        if not os.path.isdir(self.filedir):
            os.mkdir(self.filedir)
            self.download_yields()
        else:
            if not all(
                [os.path.isfile(self.filedir + os.sep + file) for file in self.files]
            ):
                self.download_yields()

        self.elements, self.atomic_num = self.get_element_list()

        # Asymptotic giant branch
        agb_yld = self.load_agb_yields()
        self.agb = Source(self.elements, [self.metal, self.mass], agb_yld)

    def agb_yields(
        self,
        elements: List[str],
        starPop: StarPopulation,
        interpolate: str = "nearest",
        extrapolate: bool = False,
    ) -> Quantity["mass"]:
        """Interpolate yields from stars on the asymptotic giant branch for specified elements.
        Stellar parameters can be provided as single value, array + single value, or arrays.

        Args:
            elements: list of elements, as specified by symbols (e.g., ['H'] for hydrogen).
            starPop: stellar population object
            interpolate: passed as method to scipy.interpolate.RegularGridInterpolator
            extrapolate: if False, then mass, metal, and rot are set to table if outside bound.
        Returns:
            List of yields matching provided element list

        """
        args = [
            starPop["metals"].value,
            starPop["mass"].to(u.M_sun).value,
        ]
        return (
            self.agb.get_yld(
                elements, args, interpolate=interpolate, extrapolate=extrapolate
            )
            * u.M_sun
        )

    def get_element_list(self) -> None:
        """Read element symbols from tables, atomic numbers not available."""

        elements = []
        atomic_numbers = []

        data = np.genfromtxt(
            "/Users/eric/Projects/arsenal-popsynth/arsenal_gear/arsenal_gear/element_yields/NuGrid/Battino2019_Yields_table.txt",
            skip_header=4,
            skip_footer=5,
            dtype=[("element", "U10")],
            usecols=[0],
            encoding="utf-8",
        )
        data["element"] = np.array(
            [re.sub(r"[^A-Za-z]", "", s) for s in data["element"]]
        )
        vals, idx = np.unique(data["element"], return_index=True)
        elements = [str(s) for s in vals[np.argsort(idx)]]

        return elements, atomic_numbers

    def load_agb_yields(self) -> np.ndarray:
        """Load tables of yields ejected during agb."""

        agb_yld = np.zeros([len(self.elements), self.metal.size, self.mass.size])

        data = np.genfromtxt(
            self.filedir + os.sep + self.files[0],
            skip_header=4,
            skip_footer=5,
            usecols=[0, 1, 2, 3, 4, 5, 6],
            dtype=[
                ("element", "U10"),
                ("m2z1m2", float),
                ("m3z1m2", float),
                ("m2z2m2", float),
                ("m3z2m2", float),
                ("m2z3m2", float),
                ("m3z3m2", float),
            ],
            encoding="utf-8",
        )
        data["element"] = np.array(
            [re.sub(r"[^A-Za-z]", "", s) for s in data["element"]]
        )

        for model in data.dtype.names[1:]:
            (ind_met, ind_mass) = self._get_index_from_model(model)
            agb_yld[:, ind_met, ind_mass] = np.atleast_1d(
                [
                    data[model][(data["element"] == element)].sum()
                    for element in self.elements
                ]
            )

        data = np.genfromtxt(
            self.filedir + os.sep + self.files[1],
            skip_header=4,
            skip_footer=5,
            usecols=[0, 1, 2],
            dtype=[("element", "U10"), ("m3z1m3", float), ("m2z1m3", float)],
            encoding="utf-8",
        )
        data["element"] = np.array(
            [re.sub(r"[^A-Za-z]", "", s) for s in data["element"]]
        )

        for model in data.dtype.names[1:]:
            (ind_met, ind_mass) = self._get_index_from_model(model)
            agb_yld[:, ind_met, ind_mass] = np.atleast_1d(
                [
                    data[model][(data["element"] == element)].sum()
                    for element in self.elements
                ]
            )

        return agb_yld

    def download_yields(self) -> None:
        """Downloader for tabulated yield files."""
        from ..utils.scraper import downloader

        print(
            "NuGrid (Battinoi et al., 2019, 2021) tables not found, starting download."
        )
        downloader(
            self.filedir + os.sep + self.files[0],
            self.b19_url + "/Yields_table.txt",
            message=None,
        )
        downloader(
            self.filedir + os.sep + self.files[1],
            self.b21_url + "/Yields_Table.txt",
            message=None,
        )

    @staticmethod
    def _get_index_from_model(model: str) -> int:
        """Convenience function for converting table rotation labels into table index."""
        ind_mods = {
            "m2z1m3": (0, 0),
            "m2z1m2": (1, 0),
            "m2z2m2": (2, 0),
            "m2z3m2": (3, 0),
            "m3z1m3": (0, 1),
            "m3z1m2": (1, 1),
            "m3z2m2": (2, 1),
            "m3z3m2": (3, 1),
        }
        if model in ind_mods:
            return ind_mods[model]
        else:
            raise ValueError("Model does not exist.")
