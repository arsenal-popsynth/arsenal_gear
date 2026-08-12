"""
data_converter.py
================================

This file defines the interface to various binary evolution models
through downloading, reorganzing and interpreting their outputs.
"""

import os
from abc import ABC, abstractmethod
from multiprocessing.pool import ThreadPool as Pool
from pathlib import Path

import numpy as np
import pandas as pd
import tqdm

from .be_data_structures import BinaryStarTrackSet, SingleStarTrackSet

os.environ["OPENBLAS_NUM_THREADS"] = "1"


class BinaryEvolutionConverter(ABC):
    """
    Abstract base class for converting binary evolution models from various
    sources to a standard Arsenal binary evolution format.
    """

    def __init__(self, **kwargs) -> None:
        # [Fe/H]
        self.met = kwargs.get("met", 0.014)
        # Directories to read and write data
        self.input_dir = kwargs.get("input_dir", None)
        self.output_dir = kwargs.get("output_dir", None)
        # Output times
        self.overwrite = kwargs.get("overwrite", False)

    @abstractmethod
    def convert_single_data(self) -> SingleStarTrackSet:
        """
        Abstract method for converting single star evolutionary tracks
        into an Arsenal stellar evolution track set.
        """

    @abstractmethod
    def convert_binary_data(self) -> BinaryStarTrackSet:
        """
        Abstract method for converting binary evolutionary tracks into
        an Arsenal binary evolution track set.
        """


class BPASSConverter(BinaryEvolutionConverter):
    """
    Class for reading BPASS data and converting it to the Arsenal binary evolution
    format. This is an instantiation of the BinaryEvolutionConverter base class.
    """

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

    def __init__(self, **kwargs) -> None:
        """
        Args:
            kwargs: Keyword arguments for the binary evolution table.

        Methods:
            convert_single_data     Processes single stellar track data into a SingleStarTrackSet
            convert_binary_data     Processes binary stellar track data into a BinaryStarTrackSet
        """
        # set input parameters
        super().__init__(**kwargs)

        if self.met >= 1e-3:
            self.metstr = "z" + str(int(self.met * 1000)).zfill(3)
        else:
            self.metstr = "zem" + str(-1 * int(np.log10(self.met)))
        if self.metstr not in self.mets:
            raise ValueError("Metallicity must be one of: " + str(self.mets))

        # Consistent format for directories
        if self.input_dir[-1] == "/":
            self.input_dir: str = self.input_dir
        else:
            self.intput_dir: str = self.input_dir + "/"
        if self.output_dir[-1] == "/":
            self.output_dir: str = self.output_dir
        else:
            self.output_dir: str = self.output_dir + "/"

    def convert_single_data(self):
        """
        Converts BPASS data for single stars into an Arsenal-readable
        SingleStarTrackSet.
        """
        # Create directory if it does not already exists
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

        model_directory = self.input_dir + "NEWSINMODS/" + self.metstr
        files = []

        # Scan directory
        with os.scandir(model_directory) as all_models:
            for model in all_models:
                if model.is_file() and model.name.startswith("sneplot"):
                    files.append(model.name)
        files.sort()

        # Function to extract the data
        def extract_data(model):

            data = np.genfromtxt(model_directory + "/" + model)
            model_split = model.split("-")
            d = {
                "model": str(int(round(float(model_split[-1]) * 100))).zfill(
                    5
                ),  # model name as 100*M
                "time": data[:, 1].astype("float"),  # time in yr
                "mass": data[:, 5].astype("float"),  # mass in MSun
                "logL": data[:, 4].astype("float"),  # log Lbol in Lsun
                "logT": data[:, 3].astype("float"),  # log Teff in K
                "logR": data[:, 2].astype("float"),  # log R in Rsun
            }

            small_df = pd.DataFrame(data=d)

            return small_df

        pool = Pool()
        results = list(tqdm.tqdm(pool.map(extract_data, files), total=len(files)))

        frames = []
        models = []

        for i in range(len(results)):
            frames.append(results[i])
            models.append((results[i].model.values[0]).zfill(5))

        pool.close()

        # Sort the data by model name with leading zeros to ensure correct order
        sorted_indices = np.argsort(models)
        frames = [frames[i] for i in sorted_indices]
        data = pd.concat(frames, ignore_index=True)

        if ("singles_" + self.metstr + ".pkl.gz") not in os.listdir(
            self.output_dir
        ) or self.overwrite:
            print("Saving processed data to", self.output_dir)

            data.to_pickle(
                self.output_dir + "/singles_" + self.metstr + ".pkl.gz",
                compression="gzip",
            )

        else:
            print("Cannot save model. Try setting overwrite=True...")

        return

    def convert_binary_data(self):
        """
        Converts BPASS data for binary stars into an Arsenal-readable
        BinaryStarTrackSet.
        """

        # Load single star data for secondary star evolution
        singles_fname = self.output_dir + "/singles_" + self.metstr + ".pkl.gz"
        if not os.path.exists(singles_fname):
            raise FileNotFoundError(
                f"Single star data file '{singles_fname}' not found. Please run convert_single_data() first."
            )
        else:
            singles = pd.read_pickle(
                self.output_dir + "/singles_" + self.metstr + ".pkl.gz",
                compression="gzip",
            )

        single_masses = np.empty(len(singles.model.values))
        for i in range(len(single_masses)):
            single_masses[i] = int(round(float(singles.model.values[i]))) / 100

        # Create directory if it does not already exists
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

        model_directory = self.input_dir + "/NEWBINMODS/NEWBINMODS/" + self.metstr

        files = []
        # Scan directory
        with os.scandir(model_directory) as all_models:
            for model in all_models:
                if model.is_file() and model.name.startswith("sneplot"):
                    files.append(model.name)
        files.sort()

        # Function to extract the data
        def extract_data(model):

            data = np.genfromtxt(model_directory + "/" + model)
            model_split = model.split("-")

            d = {
                "model": str(int(float(model_split[-3]) * 100)).zfill(5)
                + "_"
                + str(int(float(model_split[-2]) * 100)).zfill(3)
                + "_"
                + str(int(float(model_split[-1]) * 100)).zfill(
                    3
                ),  # model name as 100*M
                "time": data[:, 1].astype("float"),  # time in yr
                "p_mass": data[:, 5].astype("float"),  # mass in MSun
                "p_logL": data[:, 4].astype("float"),  # log Lbol in Lsun
                "p_logT": data[:, 3].astype("float"),  # log Teff in K
                "p_logR": data[:, 2].astype("float"),  # log R in Rsun
                "s_mass": data[:, 37].astype("float"),  # mass in MSun
                "s_logL": data[:, 48].astype("float"),  # log Lbol in Lsun
                "s_logT": data[:, 47].astype("float"),  # log Teff in K
                "s_logR": data[:, 46].astype("float"),  # log R in Rsun
            }

            combined_df = pd.DataFrame(data=d)

            p_dm = combined_df.p_mass.values[1:] - combined_df.p_mass.values[:-1]
            s_dm = combined_df.s_mass.values[1:] - combined_df.s_mass.values[:-1]
            # Check for mergers
            merged = np.where((p_dm > 0) & (s_dm == 0))[0]
            if len(merged) > 0:
                merger = merged[0]
                merged_star = {
                    "model": combined_df.model.values[merger:],
                    "time": combined_df.time.values[merger:],
                    "s_mass": np.zeros(len(combined_df.time.values[merger:])).astype(
                        "float"
                    ),
                    "s_logL": np.zeros(len(combined_df.time.values[merger:])).astype(
                        "float"
                    ),
                    "s_logT": np.zeros(len(combined_df.time.values[merger:])).astype(
                        "float"
                    ),
                    "s_logR": np.zeros(len(combined_df.time.values[merger:])).astype(
                        "float"
                    ),
                }
                merger_df = pd.DataFrame(
                    data=merged_star,
                    index=range(merger, merger + len(combined_df.time.values[merger:])),
                )
                combined_df.update(merger_df)
            else:
                # Get companion effective mass
                m_init = combined_df.s_mass.values[0]
                m_max = np.max(combined_df.s_mass.values)
                if (m_max > m_init) and (m_max >= 2):
                    m_eff = m_max
                else:
                    m_eff = m_init
                # Make sure to select this system
                t_ind = np.where(combined_df.s_mass.values == m_max)[0][0]
                t_eff = combined_df.time.values[-1] - combined_df.time.values[t_ind]
                # Match to model
                m_ind = np.argmin(np.abs(single_masses - m_eff))
                _star = np.where(
                    singles.model
                    == str(int(round(single_masses[m_ind] * 100))).zfill(5)
                )[0]
                _time = np.argmin(np.abs(singles.time[_star] - t_eff))

                evolved_star = {
                    "model": combined_df.model.values[t_ind],
                    "time": singles.time.values[_star[0] : _star[0] + _time + 1]
                    + t_eff,
                    "p_mass": combined_df.p_mass.values[-1],
                    "p_logL": float("NaN"),
                    "p_logT": float("NaN"),
                    "p_logR": float("NaN"),
                    "s_mass": singles.mass.values[
                        _star[0] : _star[0] + _time + 1
                    ].astype("float"),
                    "s_logL": singles.logL.values[
                        _star[0] : _star[0] + _time + 1
                    ].astype("float"),
                    "s_logT": singles.logT.values[
                        _star[0] : _star[0] + _time + 1
                    ].astype("float"),
                    "s_logR": singles.logR.values[
                        _star[0] : _star[0] + _time + 1
                    ].astype("float"),
                }
                evolved_df = pd.DataFrame(data=evolved_star)
                pd.concat([combined_df, evolved_df], ignore_index=True)

            return combined_df

        pool = Pool()
        results = list(tqdm.tqdm(pool.map(extract_data, files), total=len(files)))

        frames = []
        models = []

        for i in range(len(results)):
            frames.append(results[i])
            models.append((results[i].model.values[0]).zfill(5))

        pool.close()

        # Sort the data by model name with leading zeros to ensure correct order
        sorted_indices = np.argsort(models)
        frames = [frames[i] for i in sorted_indices]
        data = pd.concat(frames, ignore_index=True)

        ## Meets the criteria for rejuvenation
        # if (m_eff >= 1.05 * data[system, 5, 0]) and (m_eff > 2):
        # m_closest = m_single[np.argmin(np.abs(m_single - m_eff))]

        if ("binaries_" + self.metstr + ".pkl.gz") not in os.listdir(
            self.output_dir
        ) or self.overwrite:
            print("Saving processed data to", self.output_dir)

            data.to_pickle(
                self.output_dir + "/binaries_" + self.metstr + ".pkl.gz",
                compression="gzip",
            )

        else:
            print("Cannot save model. Try setting overwrite=True...")

        return


class MPAConverter(BinaryEvolutionConverter):
    """
    Class for reading MPA/Bonn stellar model data and converting it to the Arsenal binary
    evolution format. This is an instantiation of the BinaryEvolutionConverter base class.
    """

    mets = [
        "MW",
        "LMC",
        "SMC",
    ]

    def __init__(self, **kwargs) -> None:
        """
        Args:
            kwargs: Keyword arguments for the binary evolution table.

        Methods:
            convert_single_data     Processes single stellar track data into a SingleStarTrackSet
            convert_binary_data     Processes binary stellar track data into a BinaryStarTrackSet
        """
        # set input parameters
        super().__init__(**kwargs)

        if self.met not in self.mets:
            raise ValueError("Metallicity must be one of: " + str(self.mets))

        # Consistent format for directories
        if self.input_dir[-1] == "/":
            self.input_dir: str = self.input_dir
        else:
            self.intput_dir: str = self.input_dir + "/"
        if self.output_dir[-1] == "/":
            self.output_dir: str = self.output_dir
        else:
            self.output_dir: str = self.output_dir + "/"

    def convert_single_data(self):
        """
        Converts MPA/Bonn stellar model data for single stars into an Arsenal-readable
        SingleStarTrackSet.
        """

        # Create directory if it does not already exists
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

        model_directory = self.input_dir + "single_" + self.met
        files = []

        # Scan directory
        with os.scandir(model_directory) as all_models:
            for model in all_models:
                if model.is_file() and model.name.endswith("_0_compressed.pkl.gz"):
                    files.append(model.name)
        files.sort()

        # Function to extract the data
        def extract_data(model):

            df = pd.read_pickle(model_directory + "/" + model, compression="gzip")

            d = {
                "model": str(int(round(10 ** float(model[:5]) * 100))).zfill(
                    5
                ),  # model name as 100*M
                "time": df.star_age.values.astype("float"),  # time in yr
                "mass": df.star_mass.values.astype("float"),  # mass in MSun
                "logL": df.log_L.values.astype("float"),  # log Lbol in Lsun
                "logT": df.log_Teff.values.astype("float"),  # log Teff in K
                "logR": df.log_R.values.astype("float"),  # log R in Rsun
            }

            small_df = pd.DataFrame(data=d)

            return small_df

        pool = Pool()
        results = list(tqdm.tqdm(pool.map(extract_data, files), total=len(files)))

        frames = []

        for i in range(len(results)):
            frames.append(results[i])

        pool.close()

        data = pd.concat(frames, ignore_index=True)

        if ("singles_" + self.met + ".pkl.gz") not in os.listdir(
            self.output_dir
        ) or self.overwrite:
            print("Saving processed data to", self.output_dir)

            data.to_pickle(
                self.output_dir + "/singles_" + self.met + ".pkl.gz", compression="gzip"
            )

        else:
            print("Cannot save model. Try setting overwrite=True...")

        return

    def convert_binary_data(self):
        """
        Converts MPA/Bonn stellar model data for binary stars into an Arsenal-readable
        BinaryStarTrackSet.
        """

        # Create directory if it does not already exists
        # Path(self.output_dir).mkdir(parents=True, exist_ok=True)

        primary_directory = self.input_dir + self.met + "/primary/"
        secondary_directory = self.input_dir + self.met + "/secondary/"
        files = []

        # Scan directory
        with os.scandir(primary_directory) as subdirectories:
            for subdirectory in subdirectories:
                with os.scandir(
                    primary_directory + "/" + subdirectory.name + "/"
                ) as all_models:
                    for model in all_models:
                        if model.is_file() and model.name.endswith(
                            "_compressed.pkl.gz"
                        ):
                            files.append(subdirectory.name + "/" + model.name)
        files.sort()

        # Function to extract the data
        def extract_data(model):

            df_1 = pd.read_pickle(
                primary_directory + "/" + model,
                compression="gzip",
            )
            df_2 = pd.read_pickle(
                secondary_directory + "/" + model,
                compression="gzip",
            )

            pad = len(df_2.star_age.values) - len(df_1.star_age.values)

            if pad > 0:

                d = {
                    "model": str(int(round(10 ** float(model[6:11]) * 100))).zfill(5)
                    + "_"
                    + str(int(float(model[12:17]) * 100)).zfill(3)
                    + "_"
                    + str(int(float(model[18:23]) * 100)).zfill(3),
                    "time": df_2.star_age.values.astype("float"),  # time in yr
                    "p_mass": np.concatenate(
                        (
                            df_1.star_mass.values.astype("float"),
                            df_1.star_mass.values.astype("float")[-1] * np.ones(pad),
                        )
                    ),
                    "p_logL": np.concatenate(
                        (df_1.log_L.values.astype("float"), float("NaN") * np.ones(pad))
                    ),
                    "p_logT": np.concatenate(
                        (
                            df_1.log_Teff.values.astype("float"),
                            float("NaN") * np.ones(pad),
                        )
                    ),
                    "p_logR": np.concatenate(
                        (df_1.log_R.values.astype("float"), float("NaN") * np.ones(pad))
                    ),
                    "s_mass": df_2.star_mass.values.astype("float"),  # mass in MSun
                    "s_logL": df_2.log_L.values.astype("float"),  # log Lbol in Lsun
                    "s_logT": df_2.log_Teff.values.astype("float"),  # log Teff in K
                    "s_logR": df_2.log_R.values.astype("float"),  # log R in Rsun
                }

            # If merger
            else:

                d = {
                    "model": str(int(round(10 ** float(model[6:11]) * 100))).zfill(5)
                    + "_"
                    + str(int(float(model[12:17]) * 100)).zfill(3)
                    + "_"
                    + str(int(float(model[18:23]) * 100)).zfill(3),
                    "time": df_1.star_age.values.astype("float"),  # time in yr
                    "p_mass": df_1.star_mass.values.astype("float"),  # mass in MSun
                    "p_logL": df_1.log_L.values.astype("float"),  # log Lbol in Lsun
                    "p_logT": df_1.log_Teff.values.astype("float"),  # log Teff in K
                    "p_logR": df_1.log_R.values.astype("float"),  # log R in Rsun
                    "s_mass": np.concatenate(
                        (
                            df_2.star_mass.values.astype("float"),
                            float("NaN") * np.ones(-1 * pad),
                        )
                    ),
                    "s_logL": np.concatenate(
                        (
                            df_2.log_L.values.astype("float"),
                            float("NaN") * np.ones(-1 * pad),
                        )
                    ),
                    "s_logT": np.concatenate(
                        (
                            df_2.log_Teff.values.astype("float"),
                            float("NaN") * np.ones(-1 * pad),
                        )
                    ),
                    "s_logR": np.concatenate(
                        (
                            df_2.log_R.values.astype("float"),
                            float("NaN") * np.ones(-1 * pad),
                        )
                    ),
                }

            combined_df = pd.DataFrame(data=d)

            return combined_df

        pool = Pool()
        results = list(tqdm.tqdm(pool.map(extract_data, files), total=len(files)))

        frames = []

        for i in range(len(results)):
            frames.append(results[i])

        pool.close()

        data = pd.concat(frames, ignore_index=True)

        if ("binaries_" + self.met + ".pkl.gz") not in os.listdir(
            self.output_dir
        ) or self.overwrite:
            print("Saving processed data to", self.output_dir)

            data.to_pickle(
                self.output_dir + "/binaries_" + self.met + ".pkl.gz",
                compression="gzip",
            )

        else:
            print("Cannot save model. Try setting overwrite=True...")

        return
