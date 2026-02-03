# arsenal_gear/element_yields/sukhbold2016.py
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from scipy.interpolate import interp1d

from ..population import StarPopulation
from .source import Source
from .yieldtables import YieldTables


class Sukhbold2016(YieldTables):
    """
    Sukhbold et al. (2016) CCSN + wind yields.
    Uses pre-parsed CSVs:
      - sukhbold_ccsn_yields.csv
      - sukhbold_wind_yields.csv
      - sukhbold_explosion_energies.csv
    All located in ../../data/Sukhbold2016/
    """

    def __init__(self):
        super().__init__()
        base_dir = (
            Path(__file__).resolve().parent.parent.parent / "data" / "Sukhbold2016"
        )

        # Load ejecta (CCSN) and wind yields
        ccsn_path = base_dir / "sukhbold_ccsn_yields.csv"
        wind_path = base_dir / "sukhbold_wind_yields.csv"
        energy_path = base_dir / "sukhbold_explosion_energies.csv"

        if not ccsn_path.exists() or not wind_path.exists():
            raise FileNotFoundError(
                f"Missing yield CSVs in {base_dir}\n"
                "Run parse_sukhbold_all.py to generate them."
            )

        df_ccsn = pd.read_csv(ccsn_path).set_index("mass_Msun").fillna(0.0)
        df_wind = pd.read_csv(wind_path).set_index("mass_Msun").fillna(0.0)

        # Auto-sum isotopes to elements (e.g. o16_ejecta + o17_ejecta + o18_ejecta → O_EJECTA)
        def add_elemental_sums(df, suffix):
            elemental = {}
            for col in df.columns:
                if col == "mass_Msun" or "_" not in col:
                    continue
                # Get base isotope name (before _ejecta or _wind)
                base = col.rsplit("_", 1)[0]
                elem_raw = "".join(c for c in base if not c.isdigit())
                elem = elem_raw.title() if elem_raw else ""
                if elem not in elemental:
                    elemental[elem] = []
                elemental[elem].append(col)

            # Collect all summed Series
            summed_series = []
            for elem, cols in elemental.items():
                if cols:
                    summed = (
                        df[cols].sum(axis=1, skipna=True).rename(f"{elem}_{suffix}")
                    )
                    summed_series.append(summed)

            # One single concat if we have any
            if summed_series:
                summed_df = pd.concat(summed_series, axis=1)
                df = pd.concat([df, summed_df], axis=1)

            return df

        df_ccsn = add_elemental_sums(df_ccsn, "EJECTA")
        df_wind = add_elemental_sums(df_wind, "WIND")

        self.mass = df_ccsn.index.values
        self.elements = [c.rsplit("_", 1)[0] for c in df_ccsn.columns if "_" in c]

        # Build CCSN Source
        ccsn_grid = np.zeros((len(self.elements), 1, 1, len(self.mass)))
        for i, el in enumerate(self.elements):
            col = f"{el}_EJECTA"
            if col in df_ccsn:
                ccsn_grid[i, 0, 0, :] = df_ccsn[col].values
        self.ccsn = Source(
            self.elements, [np.array([0.0]), np.array([0.02]), self.mass], ccsn_grid
        )

        # Build wind Source
        wind_grid = np.zeros((len(self.elements), 1, 1, len(self.mass)))
        for i, el in enumerate(self.elements):
            col = f"{el}_WIND"
            if col in df_wind:
                wind_grid[i, 0, 0, :] = df_wind[col].values
        self.wind = Source(
            self.elements, [np.array([0.0]), np.array([0.02]), self.mass], wind_grid
        )

        # Load explosion energies
        self.has_energies = False
        if energy_path.exists():
            df_en = pd.read_csv(energy_path)
            self._energy_interp = interp1d(
                df_en["mass_Msun"],
                df_en["explosion_energy"],
                kind="nearest",
                fill_value=0,
                bounds_error=False,
            )
            self.has_energies = True
        else:
            warnings.warn(f"Explosion energies not found: {energy_path}")

    def get_explosion_energy(self, mass):
        if not self.has_energies:
            return np.full_like(np.atleast_1d(mass), np.nan, dtype=float)
        return self._energy_interp(np.atleast_1d(mass))

    def wind_yields(
        self,
        elements,
        starPop: StarPopulation,
        interpolate="nearest",
        extrapolate=False,
    ):
        args = [
            np.zeros(len(starPop["mass"])),
            np.full(len(starPop["mass"]), 0.02),
            starPop["mass"].to(u.Msun).value,
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
        starPop: StarPopulation,
        interpolate="nearest",
        extrapolate=False,
    ):
        if interpolate != "nearest":
            warnings.warn(
                "Sukhbold2016: forcing nearest-neighbor due to explodability islands"
            )
            interpolate = "nearest"

        args = [
            np.zeros(len(starPop["mass"])),
            np.full(len(starPop["mass"]), 0.02),
            starPop["mass"].to(u.Msun).value,
        ]
        yld_array = (
            self.ccsn.get_yld(
                elements, args, interpolate=interpolate, extrapolate=extrapolate
            )
            * u.M_sun
        )

        return {el: yld_array[i] for i, el in enumerate(elements)}
