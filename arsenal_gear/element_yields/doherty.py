"""
doherty
=======

This submodule contains all the code required to load yields from Doherty et al. (2014a,b).
"""

import re
import zipfile
from typing import List

import astropy.units as u
import numpy as np

from ..formation import SinglePop
from ..utils.citations import cite
from .source import Source
from .yieldtables import YieldTables


@cite("10.1093/mnras/stt1877")
class Doherty2014a(YieldTables):
    """
    Include yields from super-asymptotic giant branch stars from Doherty et al. (2014a) for majority of
    elements for stars with mass [6.5, 7.0, 7.5, 8.0, 8.5, 9.0] Msun, with initial metallicities
    [0.004, 0.008, 0.02].

    Yields available in suplementary material of Doherty et al. (2014a) but cannot be automatically downloaded.

    Reference: Doherty et al., 2014a, MNRAS, 437, 195

    Attributes:
        d14_url: yield source website
        mass: tabulated masses
        metal: tabulated metallicities
        filedir: directory of yield tables
        name: name of yields tables
        elements: elements available in table
        atomic_num: atomic numbers of available elements
        sagb: source object for stellar winds (super asymptotic giant branch)
    """

    models = ["VW93", "VW-M"]
    d14_url = (
        "https://academic.oup.com/mnras/article/441/1/582/980303#supplementary-data"
    )
    file = "998013_Supplementary_Data.zip"
    mass = np.array([6.5, 7.0, 7.5, 8.0, 8.5, 9.0])
    metal = np.array([0.004, 0.008, 0.02])

    def __init__(self, model: str = "VW93", **kwargs):
        """
        Initializes tabulated yields and their coordinates for the Doherty2014a yields.

        Arguments:
            model: model used for mass loss rates.
        """
        if model not in self.models:
            raise ValueError(f"Model {model} does not exist.")
        self.model = model

        super().__init__(**kwargs)
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
        if not self.available:
            return

        self.elements, self.atomic_num = self.get_element_list()

        # Super asymptotic giant branch
        sagb_yld = self.load_sagb_yields()
        self.sagb = Source(self.elements, [self.metal, self.mass], sagb_yld)

    def sagb_args(self, starPop: SinglePop) -> List:
        """Argument structure for Doherty (2014a) sagb tables"""
        args = [
            starPop.metallicity.value * starPop.ones_like,
            starPop.masses.to(u.M_sun).value,
        ]
        return args

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
        import warnings

        print("Doherty (2014a) tables not found, starting download.")
        warnings.warn(
            f"Automatic download for these yields is not available. "
            f"Please download yields from suplementary material in {self.name} and place"
            f"in {self.filedir} (follow link {self.d14_url} for download)"
        )
        self.available = False
        # downloader(self.filedir / self.file, self.d14_url, message=None)


@cite("10.1093/mnras/stu571")
class Doherty2014b(YieldTables):
    """
    Include yields from super-asymptotic giant branch stars from Doherty et al. (2014b) for majority of
    elements for stars with mass [6.5, 7.0, 7.5] Msun, with initial metallicities [0.0001, 0.001].

    Yields available in suplementary material of Doherty et al. (2014b) but cannot be automatically downloaded

    Reference: Doherty et al., 2014b, MNRAS, 441, 582

    Attributes:
        d14_url: yield source website
        mass: tabulated masses
        metal: tabulated metallicities
        filedir: directory of yield tables
        name: name of yields tables
        elements: elements available in table
        atomic_num: atomic numbers of available elements
        sagb: source object for stellar winds (super asymptotic giant branch)
    """

    models = ["VW93", "B95"]
    d14_url = (
        "https://academic.oup.com/mnras/article/441/1/582/980303#supplementary-data"
    )
    file = "stu571_Supplementary_Data.zip"
    mass = np.array([6.5, 7.0, 7.5])
    metal = np.array([0.0001, 0.001])

    def __init__(self, model: str = "VW93", **kwargs):
        """
        Initializes tabulated yields and their coordinates for the Doherty2014b yields.

        Arguments:
            model: model used for mass loss rates.
        """
        if model not in self.models:
            raise ValueError(f"Model {model} does not exist.")
        self.model = model

        super().__init__(**kwargs)
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
        if not self.available:
            return

        self.elements, self.atomic_num = self.get_element_list()

        # Super asymptotic giant branch
        sagb_yld = self.load_sagb_yields()
        self.sagb = Source(self.elements, [self.metal, self.mass], sagb_yld)

    def sagb_args(self, starPop: SinglePop) -> List:
        """Argument structure for Doherty (2014b) sagb tables"""
        args = [
            starPop.metallicity.value * starPop.ones_like,
            starPop.masses.to(u.M_sun).value,
        ]
        return args

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
        import warnings

        print("Doherty (2014b) tables not found, starting download.")
        warnings.warn(
            f"Automatic download for these yields is not available. "
            f"Please download yields from suplementary material in {self.name} and place"
            f"in {self.filedir} (follow link {self.d14_url} for download)"
        )
        self.available = False
        # downloader(self.filedir / self.file, self.d14_url, message=None)
