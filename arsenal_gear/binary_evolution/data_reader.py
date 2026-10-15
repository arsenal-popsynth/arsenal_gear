"""
data_reader.py
==========

This file defines classes for acquiring stellar and binary evolution data from
binary evolution codes and converting it into tables in a uniform format:
one row per track point, with a "model" column naming the track.
Singles tables have columns model, time, mass, logL, logT, logR; binaries
tables have the same quantities for the primary and secondary, prefixed by
"p_" and "s_". The read_* methods turn these tables into TrackSet and
BinaryTrackSet data structures.
"""

from abc import ABC, abstractmethod
from multiprocessing.pool import ThreadPool
from pathlib import Path

import astropy.constants as const
import astropy.units as u
import numpy as np
import pandas as pd
import tqdm

from ..stellar_evolution.se_data_structures import StellarTrack, TrackSet
from ..utils import downloader, untar, unzip
from ..utils.citations import cite
from .be_data_structures import BinaryTrack, BinaryTrackSet

# variable names of the track quantities in StellarTrack.qs, shared by all tracks
_TRACK_NAMES = {
    "age_name": "time",
    "mass_name": "mass",
    "lteff_name": "logT",
    "llbol_name": "logL",
    "lrad_name": "logR",
    "lgrav_name": "logg",
}
# the same names for TrackSet, which has no current-mass column name
_SET_NAMES = {k: v for k, v in _TRACK_NAMES.items() if k != "mass_name"}
# log10(G Msun / Rsun^2 / (cm/s^2)), the solar surface gravity
_LOGG_SUN = np.log10((const.G * u.Msun / u.Rsun**2).to(u.cm / u.s**2).value)


def _model_name(m1: float, q: float | None = None, logp: float | None = None) -> str:
    """
    Name of a track: "MMMMM" for single stars, "MMMMM_QQQ_PPP" for binaries,
    from 100 * the initial (primary) mass, mass ratio and log10(period/days).
    """
    name = f"{round(100 * m1):05d}"
    if q is not None:
        name += f"_{round(100 * q):03d}_{round(100 * logp):03d}"
    return name


def _star_columns(prefix: str, mass, logL, logT, logR) -> dict:
    """
    Columns for one star, in the order used by all tables: mass/Msun,
    log10(Lbol/Lsun), log10(Teff/K) and log10(R/Rsun).
    """
    values = (mass, logL, logT, logR)
    return {
        prefix + key: np.asarray(value, dtype=float)
        for key, value in zip(("mass", "logL", "logT", "logR"), values)
    }


def _pad(x: np.ndarray, n: int, fill: float = np.nan) -> np.ndarray:
    """Pad x to length n with the value fill."""
    return np.concatenate((x, np.full(n - len(x), fill)))


def _parse_model_name(name: str) -> tuple[float, float | None, float | None]:
    """
    Inverse of _model_name: the initial (primary) mass, and for binaries the
    initial mass ratio and log10(period/days); None for single stars.
    """
    values = [int(part) / 100 for part in name.split("_")]
    if len(values) == 1:
        return values[0], None, None
    return values[0], values[1], values[2]


def _nondecreasing_rows(time: np.ndarray) -> np.ndarray:
    """
    Mask dropping rows whose time is earlier than an earlier row's, which
    removes the occasional backward time step in the raw model files. Rows
    with repeated times are kept: late-phase timesteps can be shorter than
    the printed precision of the time column while the star still evolves.
    """
    keep = np.ones(len(time), dtype=bool)
    keep[1:] = time[1:] >= np.maximum.accumulate(time)[:-1]
    return keep


def _logg(mass: np.ndarray, log_r: np.ndarray) -> np.ndarray:
    """
    log10(surface gravity / (cm/s^2)) from mass/Msun and log10(R/Rsun);
    NaN where the star has no (positive) mass or radius.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        logg = _LOGG_SUN + np.log10(mass) - 2 * log_r
    return np.where(mass > 0, logg, np.nan)


class BinaryDataReader(ABC):
    """
    Abstract base class for acquiring and converting data from binary evolution
    codes. Subclasses provide the source-specific parts: getting the data, and
    finding and reading single-star and binary model files.
    """

    # valid metallicity strings, default metallicity, and data directory name
    mets: list[str] = []
    default_met = None
    source: str = ""

    def __init__(self, **kwargs) -> None:
        """
        Args:
            kwargs: Keyword arguments for the reader, including
                met: metallicity, in the native form of the source
                rootdir: directory holding the source's data
                         (default: binary_evolution/data/<source>)
                output_dir: directory for the converted tables
                            (default: <rootdir>/processed)
                force_download: download the data if no local copy is found
                convert: convert the data into tables if they don't exist yet
                         (default: True); if False, stop once the raw data
                         is available
                overwrite: rebuild the converted tables even if they exist
                verbose: decides on verbose output
        """
        self.met = kwargs.get("met", self.default_met)
        self.verbose = kwargs.get("verbose", False)
        self.force_download = kwargs.get("force_download", False)
        self.convert = kwargs.get("convert", True)
        self.overwrite = kwargs.get("overwrite", False)
        # make data storage directories if they don't exist
        rootdir = kwargs.get("rootdir")
        if rootdir is None:
            rootdir = Path(__file__).parent.absolute() / "data" / self.source
        self.rootdir = Path(rootdir)
        self.output_dir = Path(kwargs.get("output_dir") or self.rootdir / "processed")
        self.rootdir.mkdir(parents=True, exist_ok=True)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.metstr = self._metstr(self.met)
        if self.metstr not in self.mets:
            raise ValueError("Metallicity must be one of: " + str(self.mets))

        self.get_data()

    def get_data(self) -> None:
        """
        Makes sure the converted tables exist: if any is missing (or overwrite
        is True), makes sure the raw data is available, then converts it.
        The raw model files are left in place.
        """
        kinds = [
            kind
            for kind in ("singles", "binaries")
            if self.overwrite or not self._table_path(kind).exists()
        ]
        if not kinds:
            print(f"{self.source} tables already available in {self.output_dir}.")
            return
        self._get_raw_data()
        if not self.convert:
            return
        converters = {
            "singles": self.convert_single_data,
            "binaries": self.convert_binary_data,
        }
        # singles first: the BPASS binaries need them
        for kind in kinds:
            converters[kind]()

    @abstractmethod
    def _get_raw_data(self) -> None:
        """
        Makes sure the source's raw model files are available locally,
        downloading them if allowed, or raises FileNotFoundError.
        """

    @staticmethod
    @abstractmethod
    def _metstr(met) -> str:
        """Converts a metallicity into the source's metallicity string."""

    @abstractmethod
    def _single_files(self) -> list[str]:
        """Names of the single-star model files."""

    @abstractmethod
    def _binary_files(self) -> list[str]:
        """Names of the binary model files."""

    @abstractmethod
    def _extract_single(self, model: str) -> pd.DataFrame:
        """Reads one single-star model file into a singles table."""

    @abstractmethod
    def _extract_binary(self, model: str) -> pd.DataFrame:
        """Reads one binary model into a binaries table."""

    def convert_single_data(self) -> None:
        """Converts the single-star models into a singles table on disk."""
        self._convert("singles", "single-star", self._single_files(), self._extract_single)

    def convert_binary_data(self) -> None:
        """Converts the binary models into a binaries table on disk."""
        self._convert("binaries", "binary", self._binary_files(), self._extract_binary)

    def _convert(self, kind: str, label: str, files: list[str], extract) -> None:
        """Reads all model files with extract and saves them as the kind table."""
        print(
            f"Converting {len(files)} {self.source} {label} models "
            f"(metallicity {self.metstr})..."
        )
        data = self._extract_all(extract, files, desc=f"{self.source} {kind}")
        self._save(data, kind)

    @staticmethod
    def _scan(directory: Path, pattern: str) -> list[str]:
        """Sorted paths, relative to directory, of the files matching pattern."""
        files = sorted(
            str(path.relative_to(directory))
            for path in directory.glob(pattern)
            if path.is_file()
        )
        if not files:
            raise FileNotFoundError(f"No files matching {pattern} in {directory}")
        return files

    @staticmethod
    def _extract_all(extract, files: list[str], desc: str = "") -> pd.DataFrame:
        """
        Reads all files in parallel with extract, and combines the results
        into one table sorted by model (keeping the time order within each).
        desc labels the progress bar.
        """
        with ThreadPool() as pool:
            frames = list(
                tqdm.tqdm(
                    pool.imap(extract, files), total=len(files), desc=desc, unit="model"
                )
            )
        print("Combining and sorting the models...")
        data = pd.concat(frames, ignore_index=True)
        return data.sort_values("model", kind="stable", ignore_index=True)

    def _table_path(self, kind: str) -> Path:
        """Path of the converted singles or binaries table."""
        return self.output_dir / f"{kind}_{self.metstr}.pkl.gz"

    def _save(self, data: pd.DataFrame, kind: str) -> None:
        """Saves a converted table, unless it exists and overwrite is False."""
        path = self._table_path(kind)
        if path.exists() and not self.overwrite:
            print(f"{path} already exists. Set overwrite=True to replace it.")
            return
        print(
            f"Saving {kind} table ({data.model.nunique()} models, "
            f"{len(data)} rows) to {path}"
        )
        data.to_pickle(path, compression="gzip")

    def _load(self, kind: str) -> pd.DataFrame:
        """Loads a converted table."""
        path = self._table_path(kind)
        if not path.exists():
            raise FileNotFoundError(
                f"'{path}' not found. Run convert_{kind[:-1]}_data() first."
            )
        return pd.read_pickle(path, compression="gzip")

    @staticmethod
    def _models(table: pd.DataFrame):
        """
        Yields (model name, rows) for each model of a converted table, dropping
        rows that go back in time.
        """
        for name, rows in table.groupby("model", sort=False):
            yield name, rows[_nondecreasing_rows(rows.time.to_numpy())]

    @staticmethod
    def _track(rows: pd.DataFrame, prefix: str, m0: float) -> StellarTrack:
        """
        StellarTrack of one star from its rows in a converted table, where its
        columns start with prefix ("" for single stars, "p_" or "s_" in binaries).
        """
        qs = {"time": rows.time.to_numpy(float)}
        for col in ("mass", "logL", "logT", "logR"):
            qs[col] = rows[prefix + col].to_numpy(float)
        qs["logg"] = _logg(qs["mass"], qs["logR"])
        return StellarTrack(m0=m0 * u.Msun, eeps=None, elems=[], qs=qs, **_TRACK_NAMES)

    @staticmethod
    def _set_fields(tracks: list) -> dict:
        """TrackSet fields shared by single-star and binary track sets."""
        ages = [track.age.to(u.yr).value for track in tracks]
        return {
            "min_ages": np.array([age.min() for age in ages]) * u.yr,
            "max_ages": np.array([age.max() for age in ages]) * u.yr,
            "max_eep": None,
            "eep_name": None,
            "elems": [],
            **_SET_NAMES,
        }

    def read_single_data(self) -> TrackSet:
        """
        Reads the converted single-star table into a TrackSet, with one
        StellarTrack per model and no EEPs.
        """
        tracks = [
            self._track(rows, "", _parse_model_name(name)[0])
            for name, rows in self._models(self._load("singles"))
        ]
        return TrackSet(
            masses=u.Quantity([track.m0 for track in tracks]),
            hdr_list=list(tracks[0].qs),
            tracks=tracks,
            **self._set_fields(tracks),
        )

    def read_binary_data(self) -> BinaryTrackSet:
        """
        Reads the converted binaries table into a BinaryTrackSet, with one
        BinaryTrack per system. Initial parameters come from the model names;
        the initial eccentricity e0 is 0, i.e. orbits are assumed circular,
        since neither source provides it. No orbital evolution is stored yet.
        """
        binaries = []
        for name, rows in self._models(self._load("binaries")):
            m1, q, logp = _parse_model_name(name)
            binaries.append(
                BinaryTrack(
                    primary=self._track(rows, "p_", m1),
                    secondary=self._track(rows, "s_", q * m1),
                    q0=q,
                    P0=10**logp * u.d,
                    e0=0.0,
                    period_name=None,
                    eccentricity_name=None,
                )
            )
        return BinaryTrackSet(
            masses=u.Quantity([b.primary.m0 for b in binaries]),
            hdr_list=list(binaries[0].primary.qs),
            tracks=binaries,
            q0=np.array([b.q0 for b in binaries]),
            P0=u.Quantity([b.P0 for b in binaries]),
            e0=np.array([b.e0 for b in binaries]),
            **self._set_fields(binaries),
        )


@cite("10.1017/pasa.2017.51", "10.1093/mnras/sty1353")
class BPASSReader(BinaryDataReader):
    """
    Class for acquiring BPASS stellar models and converting them into tables.
    This is an instantiation of the BinaryDataReader abstract base class.
    """

    source = "BPASS"
    default_met = 0.014
    mets = [
        "zem5",
        "zem4",
        "z001",
        "z002",
        "z003",
        "z004",
        "z006",
        "z008",
        "z010",
        "z014",
        "z020",
        "z030",
        "z040",
    ]
    # Change dl=0 to dl=1 in the Dropbox link to force download
    url = "https://www.dropbox.com/scl/fo/mpuas1xh5owmdadu0vpev/h?dl=1&e=1&rlkey=7vlk7ra6kvoztzmae8wr34kmz"
    # name of the release on the server, and the only member of it that we need
    zip_name = "BPASSv2.2.1_release-07-18-Tuatara.zip"
    tar_name = "bpass-v2.2-newmodels.tar.gz"
    # directory created by extracting tar_name
    models_name = "NEWBINMODS"
    # 0-indexed data columns of mass, logL, logT and logR for each star
    primary_cols = [5, 4, 3, 2]
    secondary_cols = [37, 48, 47, 46]

    def __init__(self, **kwargs) -> None:
        """
        Args:
            kwargs: Keyword arguments for BinaryDataReader, plus
                delete_zip: delete the (~13 GB) release zip, and the models tar
                            file, once the models have been extracted from them
        """
        self.delete_zip = kwargs.get("delete_zip", True)
        # singles table and its masses, loaded by convert_binary_data
        self._singles = None
        self._single_masses = None
        super().__init__(**kwargs)

    @staticmethod
    def _metstr(met: float) -> str:
        """Converts a metallicity mass fraction Z to a BPASS string, e.g. z014."""
        if met >= 1e-3:
            return f"z{round(met * 1000):03d}"
        return f"zem{-round(np.log10(met)):01d}"

    def _get_raw_data(self) -> None:
        """
        Looks for the BPASS models and downloads/unzips/untars them if they
        are not available. The zip and tar files are deleted once extracted
        (unless delete_zip is False).
        """
        models_dir = self.rootdir / self.models_name
        tar_path = self.rootdir / self.tar_name
        zip_path = self.rootdir / self.zip_name

        if models_dir.is_dir():
            return
        if not tar_path.is_file():
            if not zip_path.is_file():
                if not self.force_download:
                    raise FileNotFoundError(
                        f"BPASS data not found in {self.rootdir}; "
                        "set force_download=True to download it"
                    )
                downloader(zip_path, self.url, "Downloading stellar models...")
            unzip(
                zip_path,
                self.rootdir,
                target_file=self.tar_name,
                delete_zip=self.delete_zip,
            )
        untar(tar_path, self.rootdir, delete_tar=self.delete_zip)
        if not models_dir.is_dir():
            raise FileNotFoundError(
                f"Extracting {tar_path} did not produce {models_dir}. "
                "The layout of the BPASS release may have changed."
            )
        print(f"BPASS data now available in {self.rootdir}.")

    def _single_dir(self) -> Path:
        return self.rootdir / "NEWSINMODS" / self.metstr

    def _binary_dir(self) -> Path:
        return self.rootdir / "NEWBINMODS" / "NEWBINMODS" / self.metstr

    def _single_files(self) -> list[str]:
        return self._scan(self._single_dir(), "sneplot*")

    def _binary_files(self) -> list[str]:
        return self._scan(self._binary_dir(), "sneplot*")

    def _extract_single(self, model: str) -> pd.DataFrame:
        data = np.genfromtxt(self._single_dir() / model)
        m1 = float(model.split("-")[-1])
        return pd.DataFrame(
            {
                "model": _model_name(m1),
                "time": data[:, 1].astype(float),  # time in yr
                **_star_columns("", *data[:, self.primary_cols].T),
            }
        )

    def convert_binary_data(self) -> None:
        """
        Converts the binary models into a binaries table on disk. Companions
        that outlive the binary model continue along single-star models, so the
        singles table (from convert_single_data) is loaded first.
        """
        print("Loading the singles table, used to continue companions...")
        self._singles = self._load("singles")
        self._single_masses = (
            np.unique(np.round(self._singles.model.values.astype(float)).astype(int))
            / 100
        )
        super().convert_binary_data()

    def _extract_binary(self, model: str) -> pd.DataFrame:
        data = np.genfromtxt(self._binary_dir() / model)
        m1, q, logp = (float(x) for x in model.split("-")[-3:])
        model_name = _model_name(m1, q, logp)

        combined_df = pd.DataFrame(
            {
                "model": model_name,
                "time": data[:, 1].astype(float),  # time in yr
                **_star_columns("p_", *data[:, self.primary_cols].T),
                **_star_columns("s_", *data[:, self.secondary_cols].T),
            }
        )

        p_dm = combined_df.p_mass.values[1:] - combined_df.p_mass.values[:-1]
        s_dm = combined_df.s_mass.values[1:] - combined_df.s_mass.values[:-1]
        # Check for mergers: the primary grows while the secondary's mass column is frozen
        merging = (p_dm > 0) & (s_dm == 0)
        if merging.any():
            first, last = np.flatnonzero(merging)[[0, -1]]
            # mass gained by the primary during merger steps, up to each row
            # (row i + 1 follows step i, hence the leading zero)
            gain = np.concatenate(([0.0], np.cumsum(p_dm * merging)))
            s_mass = combined_df.s_mass.to_numpy().copy()
            # the secondary loses what the primary gains, until it is fully absorbed
            s_mass[first : last + 1] = np.clip(
                s_mass[first : last + 1] - gain[first : last + 1], 0.0, None
            )
            s_mass[last + 1 :] = 0.0
            combined_df["s_mass"] = s_mass
            # the secondary no longer exists as a separate star
            combined_df.loc[first:, ["s_logL", "s_logT", "s_logR"]] = np.nan
            return combined_df

        # Companion's effective mass: rejuvenated if it accreted > 5% and is >= 2 Msun
        m2_init = combined_df.s_mass.values[0]
        m2_max = np.max(combined_df.s_mass.values)
        rejuvenated = (m2_max > 1.05 * m2_init) and (m2_max >= 2)
        m2_eff = m2_max if rejuvenated else m2_init
        # the companion is treated as born when it first reaches m2_eff
        i_birth = np.flatnonzero(combined_df.s_mass.values == m2_eff)[0]
        t_birth = combined_df.time.values[i_birth]
        # companion's age at the end of the binary model
        age_end = combined_df.time.values[-1] - t_birth
        # continue with the nearest-mass single-star model, strictly after that age
        masses = self._single_masses
        m_single = masses[np.argmin(np.abs(masses - m2_eff))]
        track = self._singles[self._singles.model == _model_name(m_single)]
        track = track[track.time.values > age_end]
        evolved_df = pd.DataFrame(
            {
                "model": model_name,
                # companion age -> system time, so the join moves forward in time
                "time": t_birth + track.time.values,
                "p_mass": combined_df.p_mass.values[-1],
                "p_logL": np.nan,
                "p_logT": np.nan,
                "p_logR": np.nan,
                **_star_columns(
                    "s_", track.mass, track.logL, track.logT, track.logR
                ),
            }
        )
        return pd.concat([combined_df, evolved_df], ignore_index=True)


class MPAReader(BinaryDataReader):
    """
    Class for acquiring MPA/Bonn stellar models and converting them into tables.
    This is an instantiation of the BinaryDataReader abstract base class.
    """

    source = "MPA"
    default_met = "MW"
    mets = ["MW", "LMC", "SMC"]
    url = "https://wwwmpa.mpa-garching.mpg.de/stellgrid/data/BonnGrids/"
    # model columns of mass, logL, logT and logR
    star_cols = ["star_mass", "log_L", "log_Teff", "log_R"]

    def __init__(self, **kwargs) -> None:
        """
        Args:
            kwargs: Keyword arguments for BinaryDataReader, plus
                username, password: credentials for the MPA/Bonn data server,
                                    only needed to download
        """
        self.username = kwargs.get("username")
        self.password = kwargs.get("password")
        super().__init__(**kwargs)

    @staticmethod
    def _metstr(met: str) -> str:
        """MPA/Bonn metallicities are already strings: MW, LMC or SMC."""
        return met

    def _model_sets(self) -> list[str]:
        """Model sets on the server for this metallicity."""
        return [
            f"{self.met}/primary",
            f"{self.met}/secondary",
            f"single_{self.met}/single_{self.met}",
        ]

    def _have_raw_data(self) -> bool:
        return self._single_dir().is_dir() and self._primary_dir().is_dir()

    def _get_raw_data(self) -> None:
        """
        Looks for the MPA/Bonn models and downloads the tarballs if they are not
        available. Extracting them into place is still manual (see #46).
        """
        if self._have_raw_data():
            return
        if not self.force_download:
            raise FileNotFoundError(
                f"MPA/Bonn data not found in {self.rootdir}; "
                "set force_download=True to download it"
            )
        if self.username is None:
            raise ValueError("Downloading MPA/Bonn data needs a username and password")
        for model_set in self._model_sets():
            fname = self.rootdir / f"{model_set}.tar"
            fname.parent.mkdir(parents=True, exist_ok=True)
            downloader(
                fname,
                self.url + model_set + ".tar",
                "Downloading stellar models...",
                self.username,
                self.password,
            )
        if not self._have_raw_data():
            raise FileNotFoundError(
                f"MPA/Bonn tarballs downloaded to {self.rootdir}. Extract them into "
                f"single_{self.met}/ and {self.met}/primary/, {self.met}/secondary/ "
                "before converting (see #46)."
            )

    def _single_dir(self) -> Path:
        return self.rootdir / f"single_{self.met}"

    def _primary_dir(self) -> Path:
        return self.rootdir / self.met / "primary"

    def _secondary_dir(self) -> Path:
        return self.rootdir / self.met / "secondary"

    def _single_files(self) -> list[str]:
        return self._scan(self._single_dir(), "*_0_compressed.pkl.gz")

    def _binary_files(self) -> list[str]:
        return self._scan(self._primary_dir(), "*/*_compressed.pkl.gz")

    def _extract_single(self, model: str) -> pd.DataFrame:
        df = pd.read_pickle(self._single_dir() / model, compression="gzip")
        # the file name starts with log10(M)
        return pd.DataFrame(
            {
                "model": _model_name(10 ** float(model[:5])),
                "time": df.star_age.to_numpy(float),  # time in yr
                **_star_columns("", *(df[col] for col in self.star_cols)),
            }
        )

    def _extract_binary(self, model: str) -> pd.DataFrame:
        # the file name stores log10(M1), q and logP
        model_name = _model_name(
            10 ** float(model[6:11]), float(model[12:17]), float(model[18:23])
        )
        df_1 = pd.read_pickle(self._primary_dir() / model, compression="gzip")
        df_2 = pd.read_pickle(self._secondary_dir() / model, compression="gzip")

        # both stars share the time grid of the longer model; the shorter one is
        # padded with NaN, except that the primary keeps its final mass
        n = max(len(df_1), len(df_2))
        longer = df_2 if len(df_2) > len(df_1) else df_1
        p_cols = [df_1[col].to_numpy(float) for col in self.star_cols]
        s_cols = [df_2[col].to_numpy(float) for col in self.star_cols]
        return pd.DataFrame(
            {
                "model": model_name,
                "time": longer.star_age.to_numpy(float),  # time in yr
                **_star_columns(
                    "p_",
                    _pad(p_cols[0], n, p_cols[0][-1]),
                    *(_pad(col, n) for col in p_cols[1:]),
                ),
                **_star_columns("s_", *(_pad(col, n) for col in s_cols)),
            }
        )
