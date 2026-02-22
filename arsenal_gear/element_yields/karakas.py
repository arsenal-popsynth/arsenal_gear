"""
karakas
==========

This submodule contains all the code required to load yields from
Karakas (2010) and Karakas (2016).
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


class Karakas2010(YieldTables):
    """
    Yields available at https://cdsarc.cds.unistra.fr/viz-bin/cat/J/MNRAS/403/1413#/browse

    Reference: Karakas, 2010, MNRAS, 403, 1413
    """

    k10_url = "https://cdsarc.cds.unistra.fr/viz-bin/nph-Cat/tar.gz?J/MNRAS/403/1413"
    file = "J_MNRAS_403_1413.tar.gz"
    mass = np.array(
        [1.0, 1.25, 1.5, 1.75, 1.9, 2.0, 2.25, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]
    )
    metal = np.array([0.0001, 0.004, 0.008, 0.02])

    def __init__(self):
        """
        Usage:
            >> k10 = arsenal_gear.element_yields.Karakas2010()
            >> mass = np.linspace(8, 120, 100)*u.M_sun
            >> metals = u.dimensionless_unscaled * 0.1 * np.ones(100)
            >> tform = u.Myr * np.zeros(100)
            >> stars = arsenal_gear.population.StarPopulation(mass=mass, metals=metals, tform=tform)
            >> plt.plot(mass, k10.ccsn_yields('H', stars, interpolate='nearest'), '-')

        Attributes:
            k10_url          Yield source website
            mass             Tabulated masses
            metal            Tabulated metallicities
            filedir          Directory of yield tables
            name             Name of yields tables
            elements         Elements available in table
            atomic_num       Atomic numbers of available elements
            agb              Source object for stellar winds (asymptotic giant branch)
        """

        super().__init__()
        self.filedir = self.filedir / "data" / "Karakas"
        self.name = "Karakas (2010)"

        if not self.filedir.is_dir():
            self.filedir.mkdir(parents=True, exist_ok=True)
            self.download_yields()
        else:
            if not (self.filedir / self.file).is_file():
                self.download_yields()
        self.filedir = self.filedir / self.file

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
        yld_array = self.agb.get_yld(
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

        with tarfile.open(self.filedir, "r:*") as tar:
            agbfile = tar.extractfile("./tablea2.dat")
            lines = io.TextIOWrapper(agbfile, encoding="utf-8").readlines()

        for line in lines[3:]:
            if line.split()[3] == "g":
                return elements, atomic_num
            if line.split()[3] == "d":
                element = "H"
            else:
                element = "".join(re.findall(r"[a-zA-Z]", line.split()[3]))
                element = element[:1].upper() + element[1:]
            if element not in elements:
                elements.append(element)

        raise ValueError("Could not determine list of elements from tablea2.dec")

    def load_agb_yields(self) -> np.ndarray:
        """Load tables of yields ejected as winds from asymptotic giant branch stars."""

        agb_yld = np.zeros([len(self.elements), self.metal.size, self.mass.size])

        with tarfile.open(self.filedir, "r:*") as tar:
            for member in tar.getmembers():
                if any(s in member.name for s in ("tablea1", "tablea6", "ReadMe")):
                    # These files don't contain the data we want.
                    continue

                f = tar.extractfile(member.name)
                lines = io.TextIOWrapper(f, encoding="utf-8").readlines()
                ind_metal = self._get_metal_index_from_filename(member.name)
                for line in lines:
                    if any(s in line.split()[3] for s in ("g", "n")):
                        continue

                    if float(line.split()[0]) == 6.5:
                        # Skip extra mass for 0.02
                        continue

                    # Find mass index
                    ind_mass = np.argmin(np.abs(self.mass - float(line.split()[0])))

                    # Find element index
                    if line.split()[3] == "p":
                        ind_elem = 0
                        current_elem = "H"
                    else:
                        current_elem = "".join(re.findall(r"[a-zA-Z]", line.split()[3]))
                        current_elem = current_elem[:1].upper() + current_elem[1:]

                        ind_elem = (
                            0
                            if current_elem == "D"
                            else next(
                                i
                                for i, elem in enumerate(self.elements)
                                if elem == current_elem
                            )
                        )

                    agb_yld[ind_elem, ind_metal, ind_mass] += float(line.split()[6])

        return agb_yld

    def download_yields(self) -> None:
        """Downloader for tabulated yield files."""
        from ..utils.scraper import downloader

        print("Karakas (2010) tables not found, starting download.")
        downloader(self.filedir / self.file, self.k10_url, message=None)

    # helper functions for reading tables
    @staticmethod
    def _get_metal_index_from_filename(filename: str) -> int:
        """Convenience function for converting table metal labels into table index."""
        met_ind_files = {
            "./tablea2.dat": 3,
            "./tablea3.dat": 2,
            "./tablea4.dat": 1,
            "./tablea5.dat": 0,
        }
        if filename in met_ind_files:
            return met_ind_files[filename]
        else:
            raise ValueError("Model does not exist.")


class Karakas2016(YieldTables):
    """

    Yields available at

    References:
    """

    k16_url = "https://content.cld.iop.org/journals/0004-637X/825/1/26/revision1/"
    file = "apjaa23d9_table7.tar.gz"
    mass = np.array(
        [
            1.0,
            1.25,
            1.5,
            1.75,
            2.0,
            2.25,
            2.5,
            2.75,
            3.0,
            3.25,
            3.5,
            3.75,
            4.0,
            4.25,
            4.5,
            4.75,
            5.0,
            5.5,
            6.0,
            7.0,
        ]
    )
    metal = np.array([0.007, 0.014, 0.03])

    def __init__(self, include_overshoot: bool = False):
        """
        Args:
            include_overshoot: set True if using models that mixing layer overshoot

        Usage:
            >> k16 = arsenal_gear.element_yields.Karakas2010()
            >> mass = np.linspace(8, 120, 100)*u.M_sun
            >> metals = u.dimensionless_unscaled * 0.1 * np.ones(100)
            >> tform = u.Myr * np.zeros(100)
            >> stars = arsenal_gear.population.StarPopulation(mass=mass, metals=metals, tform=tform)
            >> plt.plot(mass, k16.ccsn_yields('H', stars, interpolate='nearest'), '-')

        Attributes:
            k16_url          Yield source website
            mass             Tabulated masses
            metal            Tabulated metallicities
            filedir          Directory of yield tables
            name             Name of yields tables
            elements         Elements available in table
            atomic_num       Atomic numbers of available elements
            agb              Source object for stellar winds (asymptotic giant branch)
        """

        super().__init__()
        self.filedir = self.filedir / "data" / "Karakas"
        self.name = "Karakas (2016)"
        self.include_overshoot = include_overshoot

        if not self.filedir.is_dir():
            self.filedir.mkdir(parents=True, exist_ok=True)
            self.download_yields()
        else:
            if not (self.filedir / self.file).is_file():
                self.download_yields()
        self.filedir = self.filedir / self.file

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
        yld_array = self.agb.get_yld(
            elements, args, interpolate=interpolate, extrapolate=extrapolate
        )

        if len(elements) == 1:
            return {elements[0]: yld_array * u.M_sun}
        return {element: yld_array[i] * u.M_sun for i, element in enumerate(elements)}

    ## Internal functions for loading data
    def get_element_list(self) -> None:
        """Read element symbols and atomic numbers from tables."""

        elements = []
        atomic_num = []

        with tarfile.open(self.filedir, "r:*") as tar:
            agbfile = tar.extractfile("yield_z007.dat")
            lines = io.TextIOWrapper(agbfile, encoding="utf-8").readlines()

        elements.append("H")
        atomic_num.append(1)
        for line in lines[4:]:
            if line[0] == "#":
                return elements, np.array(atomic_num)

            element = "".join(re.findall(r"[a-zA-Z]", line.split()[0]))
            element = element[:1].upper() + element[1:]
            if element not in elements:
                elements.append(element)
                atomic_num.append(int(line.split()[1]))

        raise ValueError("Could not determine list of elements from tablea2.dec")

    def load_agb_yields(self) -> np.ndarray:
        """Load tables of yields ejected as winds from asymptotic giant branch stars."""

        agb_yld = np.zeros([len(self.elements), self.metal.size, self.mass.size])

        with tarfile.open(self.filedir, "r:*") as tar:
            for member in tar.getmembers():
                if any(
                    s in member.name
                    for s in (
                        "._ReadMe.t7",
                        "._yield_z007.dat",
                        "._yield_z014.dat",
                        "._yield_z03.dat",
                        "ReadMe.t7",
                    )
                ):
                    # These files don't contain the data we want.
                    continue

                f = tar.extractfile(member.name)
                lines = io.TextIOWrapper(f, encoding="utf-8").readlines()
                ind_metal = self._get_metal_index_from_filename(member.name)
                ind_mass = None
                self.mmix = np.zeros([self.metal.size, self.mass.size])
                for line in lines:
                    if "# Initial mass =" in line:
                        if self.include_overshoot and "N_ov = 0.0" in line:
                            continue
                        if not self.include_overshoot and "N_ov = 1.0" in line:
                            continue

                        mass = float(
                            line.split("# Initial mass =")[1].split(", Z =")[0]
                        )
                        ind_mass = np.argmin(np.abs(self.mass - mass))

                        if mass > 7.0:
                            # Only in some models
                            continue
                        # Save in case someone is interested.
                        try:
                            self.mmix[ind_metal, ind_mass] = float(
                                line.split("M_mix =")[-1]
                            )
                        except ValueError:
                            self.mmix[ind_metal, ind_mass] = float(
                                line.split("M_mix =")[-1].split(", N_ov =")[0]
                            )

                    if mass > 7.0:
                        # Only in some models
                        continue

                    # Skip these
                    if line[0] == "#":
                        continue
                    if len(line.split()) == 5:
                        continue

                    ind_elem = np.argmin(np.abs(self.atomic_num - int(line.split()[1])))
                    agb_yld[ind_elem, ind_metal, ind_mass] = float(line.split()[6])

        return agb_yld

    def download_yields(self) -> None:
        """Downloader for tabulated yield files."""
        from ..utils.scraper import downloader

        print("Karakas (2016) tables not found, starting download.")
        downloader(self.filedir / self.file, self.k16_url + self.file, message=None)

    # helper functions for reading tables
    @staticmethod
    def _get_metal_index_from_filename(filename: str) -> int:
        """Convenience function for converting table metal labels into table index."""
        met_ind_files = {"yield_z03.dat": 2, "yield_z014.dat": 1, "yield_z007.dat": 0}
        if filename in met_ind_files:
            return met_ind_files[filename]
        else:
            raise ValueError("Model does not exist.")
