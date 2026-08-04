"""
sukhbold2016
==========

This submodule contains all the code necessary to load yields from Sukhbold et al. (2016)
"""

import os
import re
import tarfile
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.units import Quantity
from scipy.interpolate import interp1d

from ..formation import SinglePop
from .source import Source
from .yieldtables import YieldTables


class Sukhbold2016(YieldTables):
    """
    Sukhbold et al. (2016) CCSN + wind yields for masses varying from 9 to 120 Msun.

    Upon first initializtion of this class, the yield tables are
    downloaded. If this fails, the user may download the .tar.gz files

    Yields available at: https://wwwmpa.mpa-garching.mpg.de/ccsnarchive/data/SEWBJ_2015/

    Reference: Sukhbold T, Ertl T, Woosley SE, Brown JM & Janka HT, 2016, ApJ, 821, 38

    """

    base_url = "https://wwwmpa.mpa-garching.mpg.de/ccsnarchive/data/SEWBJ_2015/data/"
    nucleosynthesis_file = "nucleosynthesis_yields.tar.gz"
    energy_file = "explosion_results_PHOTB.tar.gz"

    def __init__(self, engine="N20"):

        super().__init__()

        # Valid engines per Sukhbold et al. (2016)
        # only these engines have nucleosynthesis yields
        # Z9.6 covers low mass (low high-mass (?)) of 9 - 12 Msun
        # and is automatically used alongside
        self.valid_engines = ["W18", "N20"]
        if engine not in self.valid_engines:
            raise ValueError(f"Engine must be one of {self.valid_engines}")
        self.engine = engine

        self.filedir = Path(__file__).resolve().parent / "Sukhbold2016"
        self.name = "Sukhbold et al. (2016)"

        if not self.filedir.exists():
            os.mkdir(self.filedir)

        for table in [self.nucleosynthesis_file, self.energy_file]:
            if not (self.filedir / table).exists():
                self.download_sukhbold_table(table)

        # Parse nucleosynthesis yields and explosion energies
        self.load_from_archives()

    def download_sukhbold_table(self, filename):
        from ..utils.scraper import downloader

        downloader(
            str(self.filedir / filename),
            f"{self.base_url}/{filename}",
            message=f"Downloading Sukhbold 2016 table: {filename}...",
        )

    def load_from_archives(self):
        """
        Parses nucleosynthesis yields (one .yield_table file for each mass and engine)
        and energies (single file per engine containing all ZAMS masses)
        """
        all_ejecta = []
        all_wind = []

        with tarfile.open(self.filedir / self.nucleosynthesis_file, "r:gz") as tar:
            for member in tar.getmembers():
                if not member.name.endswith(".yield_table"):
                    continue
                is_base = "Z9.6" in member.name
                is_engine = self.engine in member.name
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

        self.masses = df_ccsn["mass_Msun"].values

        all_cols = [c for c in df_ccsn.columns if c != "mass_Msun"]
        self.elements = [el for el in all_cols if el[0].isupper()]

        for col in self.elements:
            if col not in df_wind.columns:
                df_wind[col] = 0.0

        ccsn_grid = np.zeros((len(self.elements), 1, 1, len(self.masses)))
        wind_grid = np.zeros((len(self.elements), 1, 1, len(self.masses)))

        for i, el in enumerate(self.elements):
            ccsn_grid[i, 0, 0, :] = df_ccsn[el].values
            wind_grid[i, 0, 0, :] = df_wind[el].values

        # Sukhbold et al. (2016) is solar metallicity (0.02) and non-rotating (0)
        self.ccsn = Source(self.elements, [[0.0], [0.02], self.masses], ccsn_grid)
        self.wind = Source(self.elements, [[0.0], [0.02], self.masses], wind_grid)

        # This archive contains: explosion_results_PHOTB/results_{engine}
        all_energies = []
        with tarfile.open(self.filedir / self.energy_file, "r:gz") as tar:
            for member in tar.getmembers():
                if member.isfile() and (
                    "results_Z9.6" in member.name
                    or f"results_{self.engine}" in member.name
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

    def get_explosion_energy(self, mass):
        if not self.has_energies:
            return np.full_like(np.atleast_1d(mass), np.nan, dtype=float)
        return self._energy_interp(np.atleast_1d(mass))

    def wind_yields(
        self,
        elements,
        starPop: SinglePop,
        interpolate="nearest",
        extrapolate=False,
    ) -> dict[str, Quantity["mass"]]:
        args = [
            np.zeros(len(starPop.masses)),
            np.full(len(starPop.masses), 0.02),
            starPop.masses.to(u.Msun).value,
        ]
        yld_array = (
            self.wind.get_yld(
                elements, args, interpolate=interpolate, extrapolate=extrapolate
            )
            * u.M_sun
        )
        return {el: yld_array[i] for i, el in enumerate(elements)}

    def ccsn_yields(
        self,
        elements,
        starPop: SinglePop,
        interpolate="nearest",
        extrapolate=False,
    ) -> dict[str, Quantity["mass"]]:
        if interpolate != "nearest":
            warnings.warn(
                "Sukhbold2016: forcing nearest-neighbor due to explodability islands"
            )
            interpolate = "nearest"

        args = [
            np.zeros(len(starPop.masses)),
            np.full(len(starPop.masses), 0.02),
            starPop.masses.to(u.Msun).value,
        ]
        yld_array = (
            self.ccsn.get_yld(
                elements, args, interpolate=interpolate, extrapolate=extrapolate
            )
            * u.M_sun
        )

        return {el: yld_array[i] for i, el in enumerate(elements)}
