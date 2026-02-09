"""
sukhbold2016
==========

This submodule contains all the code necessary to load yields from Sukhbold et al. (2016)
"""

import os
import re
import tarfile
import warnings

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.units import Quantity
from scipy.interpolate import interp1d

from ..population import StarPopulation
from .source import Source
from .yieldtables import YieldTables


class Sukhbold2016(YieldTables):
    """
    Sukhbold et al. (2016) CCSN + wind yields for masses varying from 9 to 120 Msun
    at metallicity Z=0.02. Includes two different models for explosion engine, see
    Sukhbold et al. (2016).
    Z9.6 covers low mass (low high-mass (?)) of 9 - 12 Msun and is automatically
    used alongside.

    Upon first initializtion of this class, the yield tables are
    downloaded. If this fails, the user may download the .tar.gz files

    Yields available at: https://wwwmpa.mpa-garching.mpg.de/ccsnarchive/data/SEWBJ_2015/

    Reference: Sukhbold T, Ertl T, Woosley SE, Brown JM & Janka HT, 2016, ApJ, 821, 38

    """

    s16_url = "https://wwwmpa.mpa-garching.mpg.de/ccsnarchive/data/SEWBJ_2015/data/"
    models = ["W18", "N20"]
    files = ["nucleosynthesis_yields.tar.gz", "explosion_results_PHOTB.tar.gz"]
    mass = None  # Loaded later
    metal = 0.02

    def __init__(self, model="N20"):
        """
        Args:
            model: choise of model to load, see Sukhbold et al. (2016) for details

        Usage:
            >> s16 = arsenal_gear.element_yields.Sukhbold()
            >> mass = np.linspace(8, 120, 100)*u.M_sun
            >> metals = u.dimensionless_unscaled * 0.02
            >> tform = u.Myr * np.zeros(100)
            >> stars = arsenal_gear.population.StarPopulation(mass=mass, metals=metals, tform=tform)
            >> plt.plot(mass, s16.ccsn_yields('H', stars, interpolate='nearest'), '-')

        Attributes:
            s16_url          Yield source website
            models           Available models (explosion engine)
            mass             Tabulated masses
            metal            Available metallicity (0.02)
            filedir          Directory of yield tables
            files            Files that are downloaded
            name             Name of yields tables
            model            Engine used to drive explosion
            elements         Elements available in table
            atomic_num       Not available in this model
            wind             Source object for stellar winds (massive stars)
            ccsn             Source object for core-collapse SNe
            df_energy        Distribution of energies
            _energy_interp   Interpolator for energies
            has_energies     Flag for energies
        """
        if model not in self.models:
            raise ValueError(f"Model {model} does not exist.")

        super().__init__()
        self.filedir = self.filedir / "data" / "NuGrid" / "Sukhbold2016"
        self.name = "Sukhbold et al. (2016)"
        self.model = model

        if not self.filedir.exists():
            os.mkdir(self.filedir)

        for table in self.files:
            if not (self.filedir / table).exists():
                self.download_sukhbold_table(table)

        self.ccsn = None
        self.wind = None
        # Sets nucleosynthesis yields and explosion energies
        self.load_from_archives()

    def wind_yields(
        self,
        elements,
        starPop: StarPopulation,
        interpolate="nearest",
        extrapolate=False,
    ) -> dict[str, Quantity["mass"]]:
        """Interpolate yields from massive stars ejected as winds for specified elements.
        Stellar parameters can be provided as single value, array + single value, or arrays.

        Args:
            elements: list of elements, as specified by symbols (e.g., ['H'] for hydrogen)
            starPop: StarPopulation object
            interpolate: passed as method to scipy.interpolate.RegularGridInterpolator
            extrapolate: if False, then mass, metal, and rot are set to table if outside bound
        Returns:
            Dictionary of yields matching provided element list
        """
        if any(starPop["metals"] != 0.02):
            raise ValueError(f"{self.name} only include yields at Z = 0.02")

        elements = np.atleast_1d(elements)
        args = [
            starPop["mass"].to(u.Msun).value,
        ]
        yld_array = self.wind.get_yld(
            elements, args, interpolate=interpolate, extrapolate=extrapolate
        )

        if len(elements) == 1:
            return {elements[0]: yld_array * u.M_sun}
        return {element: yld_array[i] * u.M_sun for i, element in enumerate(elements)}

    def ccsn_yields(
        self,
        elements,
        starPop: StarPopulation,
        interpolate="nearest",
        extrapolate=False,
    ) -> dict[str, Quantity["mass"]]:
        """Interpolate yields from core-collapse supernovae for specified elements.
        Stellar parameters can be provided as single value, array + single value, or arrays.

        Args:
            elements: list of elements, as specified by symbols (e.g., ['H'] for hydrogen)
            starPop: StarPopulation object
            interpolate: passed as method to scipy.interpolate.RegularGridInterpolator
            extrapolate: if False, then mass, metal, and rot are set to limits if outside bound
        Returns:
            Dictionary of yields matching provided element list

        """
        if interpolate != "nearest":
            warnings.warn(
                "Sukhbold2016: forcing nearest-neighbor due to explodability islands"
            )
            interpolate = "nearest"

        if any(np.atleast_1d(starPop["metals"] != 0.02)):
            raise ValueError(f"{self.name} only include yields at Z = 0.02")

        elements = np.atleast_1d(elements)
        args = [
            starPop["mass"].to(u.Msun).value,
        ]
        yld_array = self.ccsn.get_yld(
            elements, args, interpolate=interpolate, extrapolate=extrapolate
        )

        if len(elements) == 1:
            return {elements[0]: yld_array * u.M_sun}
        return {element: yld_array[i] * u.M_sun for i, element in enumerate(elements)}

    def get_explosion_energy(self, mass):
        """Function that returns explosion energy interpolated from mass
        Args:
            mass: Array/list of masses

        Returns:
            energy: Explosion energy
        """
        if not self.has_energies:
            return np.full_like(np.atleast_1d(mass), np.nan, dtype=float)
        return self._energy_interp(np.atleast_1d(mass))

    def download_sukhbold_table(self, filename):
        """Downloader for tabulated yield files."""
        from ..utils.scraper import downloader

        downloader(
            str(self.filedir / filename),
            f"{self.s16_url}/{filename}",
            message=f"Downloading Sukhbold 2016 table: {filename}...",
        )

    def load_from_archives(self):
        """
        Parses nucleosynthesis yields (one .yield_table file for each mass and engine)
        and energies (single file per engine containing all ZAMS masses)
        """
        all_ejecta = []
        all_wind = []

        with tarfile.open(self.filedir / self.files[0], "r:gz") as tar:
            for member in tar.getmembers():
                if not member.name.endswith(".yield_table"):
                    continue
                is_base = "Z9.6" in member.name
                is_engine = self.model in member.name
                is_implosion = "implosions" in member.name

                if is_base or is_engine or is_implosion:
                    f = tar.extractfile(member)
                    mass, ej_dict, wi_dict = self._parse_yield_table_stream(
                        f, member.name
                    )

                    if mass is not None:
                        ej_dict["mass_Msun"] = mass
                        wi_dict["mass_Msun"] = mass

                        all_ejecta.append(ej_dict)
                        all_wind.append(wi_dict)

        df_ccsn = pd.DataFrame(all_ejecta).fillna(0.0)
        df_wind = pd.DataFrame(all_wind).fillna(0.0)

        df_ccsn = df_ccsn.groupby("mass_Msun").mean().reset_index()
        df_wind = df_wind.groupby("mass_Msun").mean().reset_index()

        df_ccsn = self.add_elemental_sums(df_ccsn)
        df_wind = self.add_elemental_sums(df_wind)

        df_ccsn = df_ccsn.sort_values("mass_Msun")
        df_wind = df_wind.sort_values("mass_Msun")

        self.mass = df_ccsn["mass_Msun"].values

        all_cols = [c for c in df_ccsn.columns if c != "mass_Msun"]
        self.elements = [el for el in all_cols if el[0].isupper()]

        for col in self.elements:
            if col not in df_wind.columns:
                df_wind[col] = 0.0

        ccsn_grid = np.zeros((len(self.elements), len(self.mass)))
        wind_grid = np.zeros((len(self.elements), len(self.mass)))

        for i, el in enumerate(self.elements):
            ccsn_grid[i, :] = df_ccsn[el].values
            wind_grid[i, :] = df_wind[el].values

        # Sukhbold et al. (2016) is solar metallicity (0.02) and non-rotating (0)
        self.ccsn = Source(self.elements, [self.mass], ccsn_grid)
        self.wind = Source(self.elements, [self.mass], wind_grid)

        # This archive contains: explosion_results_PHOTB/results_{engine}
        all_energies = []
        with tarfile.open(self.filedir / self.files[1], "r:gz") as tar:
            for member in tar.getmembers():
                if member.isfile() and (
                    "results_Z9.6" in member.name
                    or f"results_{self.model}" in member.name
                ):
                    f = tar.extractfile(member)
                    if f is not None:
                        en_df = self._parse_energy_stream(f)
                        all_energies.append(en_df)

        self.df_energy = pd.concat(all_energies).sort_values("mass_Msun")
        self._energy_interp = interp1d(
            self.df_energy["mass_Msun"].values,
            self.df_energy["explosion_energy"].values,
            kind="nearest",
            bounds_error=False,
            fill_value=0.0,
        )
        self.has_energies = True

    def _parse_yield_table_stream(self, f, filename):
        """Helper to process individual .yield_table file contents"""
        text = f.read().decode("utf-8")
        # Regex to find mass in strings like 's12.25.yield_table'
        mass_match = re.search(r"[s](\d+\.?\d*)", filename)
        if not mass_match:
            return None, None, None

        mass = float(mass_match.group(1))
        ej, wi = {}, {}
        for line in text.splitlines():
            parts = re.split(r"\s+", line.strip())
            ## header line starts with '[': [isotope]
            if len(parts) < 2 or parts[0].startswith("["):
                continue
            try:
                iso = parts[0].lower()
                if len(parts) == 2:
                    # implosion file: [isotope] [wind]
                    ej[iso] = 0.0
                    wi[iso] = float(parts[1])
                elif len(parts) == 3:
                    # Standard yield file: [isotope] [ejecta] [wind]
                    ej[iso] = float(parts[1])
                    wi[iso] = float(parts[2])
            except ValueError:
                continue
        return mass, ej, wi

    def _parse_energy_stream(self, f):
        """
        Parses the results_{engine} text streams from the tarball.
        Returns a DataFrame with columns ['mass_Msun', 'explosion_energy'].
        """
        # Decode the byte stream from the tarfile
        text = f.read().decode("utf-8")

        data = []
        # and capture the blocks between them
        model_matches = list(
            re.finditer(
                r"^([a-zA-Z]?\d+\.?\d*)\s*(?:\([^\)]*\))?\s*-*$", text, re.MULTILINE
            )
        )

        for i, match in enumerate(model_matches):
            # Determine where this model's text block ends
            start = match.start()
            end = (
                model_matches[i + 1].start()
                if i + 1 < len(model_matches)
                else len(text)
            )
            block = text[start:end]

            # 1. Extract Mass from the header
            mass_match = re.search(r"\d+\.?\d*", match.group(1))
            if not mass_match:
                continue
            mass = float(mass_match.group())

            # 2. Extract Explosion Energy (E_exp) in units of 'foe' (10^51 erg)
            # Some stars collapse to BHs and have 0 or no E_exp listed
            e_match = re.search(r"E_exp\s*=\s*([\d\.]+)\s*foe", block, re.IGNORECASE)

            if e_match:
                energy = float(e_match.group(1))
            else:
                # If E_exp isn't found, the star likely collapsed to a Black Hole
                energy = 0.0

            data.append({"mass_Msun": mass, "explosion_energy": energy})

        return pd.DataFrame(data)

    def add_elemental_sums(self, df):
        """ """
        elemental = {}
        for col in df.columns:
            if col == "mass_Msun" or "_" in col:
                continue

            # Extract element symbol: "o16" -> "O"
            base = "".join(c for c in col if not c.isdigit())
            elem = base.title()

            if elem not in elemental:
                elemental[elem] = []
            elemental[elem].append(col)

        for elem, cols in elemental.items():
            if cols:
                df[elem] = df[cols].sum(axis=1, skipna=True)

        return df
