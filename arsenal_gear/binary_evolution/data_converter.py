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
import xarray as xr

from .be_data_structures import BinaryStarTable, SingleStarTable


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
    def convert_single_data(self) -> SingleStarTable:
        """
        Abstract method for converting single star evolutionary tracks
        into an Arsenal stellar evolution table.
        """

    @abstractmethod
    def convert_binary_data(self) -> BinaryStarTable:
        """
        Abstract method for converting binary evolutionary tracks into
        an Arsenal binary evolution table.
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
            convert_single_data     Processes single stellar track data into a SingleStarTable
            convert_binary_data     Processes binary stellar track data into a BinaryStarTable
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
        SingleStarTable.
        """

        # Create directory if it does not already exists
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

        model_directory = self.input_dir + "/NEWSINMODS/" + self.metstr

        # Create the zero array
        # Number of files depends on the choice of mass limits for the IMF
        files = os.listdir(model_directory)
        num_files = len(files)
        # Times for time array
        num_times = int(1e4)  # Assume there are at most 1e4 outputs
        # Save properties as a function of time: (1) time, (2) mass,
        # (3) bolometric luminosity, (4) surface temperature, (5) radius
        data = np.zeros((num_files, 5, num_times))

        # Function to extract the data
        def extract_data(model):

            if model.startswith("sneplot"):

                _data = np.genfromtxt(model_directory + "/" + model)
                # Replace NaNs by 0s --> What about files without a companion?
                _data = np.nan_to_num(_data)
                # We want to extract the following values
                _num = np.arange(len(_data[:, 1]))

                # Indices for time in yr, mass in Msun, log Lbol in Lsun,
                # log Teff in K, log R in Rsun
                return_ids = [1, 5, 4, 3, 2]

                return np.transpose(_data[:, return_ids])  # log R in Rsun

        pool = Pool()
        i = 0
        for result in pool.imap_unordered(extract_data, files):
            data[i, :, : len(result[0, :])] = result
            i += 1

        pool.close()

        # Remove superfluous model numbers
        first_empty = 0
        j = 0
        while j < num_times:
            if len(np.nonzero(data[:, 1, j])[0]) > 0:
                first_empty += 1
                j += 1
            else:
                first_empty += 1
                j = num_times

        data = data[:, :, :first_empty]

        # If fewer models, fill with sensible values
        for i in range(len(data[:, 0, 0])):
            _zeros = np.where(data[i, 0, :] == 0)[0]
            if len(_zeros) > 1:  # if not longest
                data[i, 0, _zeros[1:]] = (
                    np.ones(len(_zeros) - 1) * data[i, 0, _zeros[1] - 1]
                )
                data[i, 1, _zeros[1:]] = (
                    np.ones(len(_zeros) - 1) * data[i, 1, _zeros[1] - 1]
                )
                data[i, 2, _zeros[1:]] = (
                    np.ones(len(_zeros) - 1) * data[i, 2, _zeros[1] - 1]
                )
                data[i, 3, _zeros[1:]] = (
                    np.ones(len(_zeros) - 1) * data[i, 3, _zeros[1] - 1]
                )
                data[i, 4, _zeros[1:]] = (
                    np.ones(len(_zeros) - 1) * data[i, 4, _zeros[1] - 1]
                )

        # Sort by mass
        _sort = np.argsort(data[:, 1, 0])
        data = data[_sort, :, :]

        if ("singles_" + self.metstr + ".h5") not in os.listdir(
            self.output_dir
        ) or self.overwrite:
            print("Saving processed data to", self.output_dir)

            masses = data[:, 1, 0]

            ds = xr.Dataset(
                data_vars=dict(
                    time=(["model", "step"], data[:, 0, :]),
                    mass=(["model", "step"], data[:, 1, :]),
                    log_Lbol=(["model", "step"], data[:, 2, :]),
                    log_Teff=(["model", "step"], data[:, 3, :]),
                    log_R=(["model", "step"], data[:, 4, :]),
                ),
                coords=dict(
                    model=("model", masses), step=("step", np.arange(first_empty))
                ),
                attrs=dict(
                    description="BPASS evolution data for single stars at Z="
                    + str(self.met)
                ),
            )
            ds.to_netcdf(self.output_dir + "/singles_" + self.metstr + ".h5")

        else:
            print("Cannot save model. Try setting overwrite=True...")

        return

    def convert_binary_data(self):
        """
        Converts BPASS data for binary stars into an Arsenal-readable
        BinaryStarTable.
        """

        # Create directory if it does not already exists
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

        model_directory = self.input_dir + "/NEWBINMODS/NEWBINMODS/" + self.metstr

        # Values for disrupted or modified systems, when looking for new systems
        M_single = np.concatenate(
            (np.arange(0.1, 10.0, 0.1), np.arange(10, 100, 1), np.arange(100, 325, 25))
        )

        # Create the zero array
        # Number of files depends on the choice of mass limits for the IMF
        files = os.listdir(model_directory)
        num_files = len(files)
        # Times for time array
        num_times = int(1e4)  # Assume there are at most 1e4 outputs
        num_times = int(1e4)  # Assume there are at most 1e4 outputs
        # Save properties as a function of time: (1) time,
        # (2) primary mass, (3) primary bolometric luminosity,
        # (4) primary surface temperature, (5) primary radius, (6)-(9) for the companion
        data = np.zeros((num_files, 9, num_times))

        model_orbits = np.zeros(num_files)

        # Function to extract the data
        def extract_data(model):

            if model.startswith("sneplot"):

                _data = np.genfromtxt(model_directory + "/" + model)
                # Replace NaNs by 0s --> What about files without a companion?
                _data = np.nan_to_num(_data)
                # We want to extract the following values
                _num = np.arange(len(_data[:, 1]))

                # Indices for time in yr, mass in Msun, log Lbol in Lsun,
                # log Teff in K, log R in Rsun, then same for companion
                return_ids = [1, 5, 4, 3, 2, 37, 48, 47, 46]

                model_split = model.split("-")

                return np.transpose(_data[:, return_ids]), model_split  # log R in Rsun

        pool = Pool()

        i = 0

        for result in pool.imap_unordered(extract_data, files):
            if result is not None:
                extracted_data, model_split = result
                data[i, :, : len(extracted_data[0, :])] = extracted_data
                model_orbits[i] = int(
                    float(model_split[2]) * 1e4
                    + float(model_split[3]) * 1e3
                    + float(model_split[4]) * 10
                )
                i += 1

        pool.close()

        # Sort by model
        sorted_models = np.argsort(model_orbits)[::-1]
        data = data[sorted_models, :, :]
        model_orbits = model_orbits[sorted_models]

        print("Pre-SN evolution done. Now checking for rejuvenation...")

        # Check for mergers and rejuvenation
        dM1 = data[:, 1, 1:] - data[:, 1, :-1]
        dM2 = data[:, 5, 1:] - data[:, 5, :-1]
        # Look for first instance of primary accreting without mass loss in companion
        # This is done because the companion properties are "frozen" at the start of the merger
        merged_systems, merged_steps = np.where((dM1 > 0) & (dM2 == 0))
        unique_systems, unique_ids = np.unique(merged_systems, return_index=True)
        unique_steps = merged_steps[unique_ids]
        for j in range(len(unique_systems)):
            data[unique_systems[j], 5, unique_steps[j] :] = 0
            data[unique_systems[j], 6, unique_steps[j] :] = -1 * np.inf
            data[unique_systems[j], 7, unique_steps[j] :] = -1 * np.inf
            data[unique_systems[j], 8, unique_steps[j] :] = -1 * np.inf

        # Evolution of companion after SN
        SN_systems = np.where(data[:, 1, 0] >= 8)[0]
        # Look for rejuventated systems
        # Here, use the last value for the subsequent evolution
        rejuv_systems, rejuv_steps = np.where(dM2 > 0)
        for system in SN_systems:
            # Has the system been rejuvenated?
            if system in rejuv_systems:

                step = rejuv_steps[np.where(rejuv_systems == system)][-1]
                time = data[system, 0, :][step]
                M_eff = data[system, 5, :][step]

                # Meets the criteria for rejuvenation
                if (M_eff >= 1.05 * data[system, 5, 0]) and (M_eff > 2):
                    M_closest = M_single[np.argmin(np.abs(M_single - M_eff))]

                # Does not meet the criteria; update time to match SN time
                else:
                    M_eff = data[system, 5, 0]
                    M_closest = M_single[np.argmin(np.abs(M_single - M_eff))]
                    step = np.where(data[system, 0, :] != 0)[0][-1]
                    time = data[system, 0, :][step]

            else:
                step = np.where(data[system, 0, :] != 0)[0][-1]
                time = data[system, 0, :][step]
                M_eff = data[system, 5, 0]
                if M_eff > 0.8:  # If actual companion, and not empty array
                    M_closest = M_single[np.argmin(np.abs(M_single - M_eff))]
                else:
                    M_closest = 0

            # Now use mass and time for subsequent evolution
            # Name-match the single BPASS files
            if M_closest == 0:
                pass
            else:
                if np.round(M_closest, 1) == np.round(M_closest, 0):
                    M_closest = str(int(M_closest))
                elif M_closest < 10:
                    M_closest = str(np.round(M_closest, 1))
                else:
                    M_closest = str(int(M_closest))
                c_file = (
                    self.input_dir
                    + "/NEWSINMODS/"
                    + self.metstr
                    + "/sneplot-"
                    + self.metstr
                    + "-"
                    + M_closest
                )
                c_data = np.genfromtxt(c_file)
                c_data = np.nan_to_num(c_data)
                # Number of data points
                nearest_num = np.argsort(np.abs(c_data[:, 1] - time))
                if c_data[nearest_num[0], 1] < (time - 1e5):
                    start_num = nearest_num[1]
                elif (time - 1e5) < c_data[nearest_num[0], 1] < time:
                    start_num = nearest_num[0]
                    c_data[nearest_num[0], 1] = time
                else:
                    start_num = nearest_num[0]
                _num = len(c_data[start_num:, 1])
                data[system, 0, step : step + _num] = c_data[
                    start_num:, 1
                ]  # time in yr
                data[system, 5, step : step + _num] = c_data[
                    start_num:, 37
                ]  # companion mass in Msun
                data[system, 6, step : step + _num] = c_data[
                    start_num:, 48
                ]  # companion log Lbol in Lsun
                data[system, 7, step : step + _num] = c_data[
                    start_num:, 47
                ]  # companion log Teff in K
                data[system, 8, step : step + _num] = c_data[
                    start_num:, 46
                ]  # companion log R in Rsun

        # Remove superfluous model numbers
        first_empty = 0
        j = 0
        while j < num_times:
            if len(np.nonzero(data[:, 1, j])[0]) > 0:
                first_empty += 1
                j += 1
            else:
                first_empty += 1
                j = num_times

        data = data[:, :, :first_empty]

        # If fewer models, fill with sensible values
        for i in range(len(data[:, 0, 0])):
            _zeros = np.where(data[i, 0, :] == 0)[0]
            if len(_zeros) > 1:  # if not longest
                data[i, 0, _zeros[1:]] = (
                    np.ones(len(_zeros) - 1) * data[i, 0, _zeros[1] - 1]
                )
                data[i, 1, _zeros[1:]] = (
                    np.ones(len(_zeros) - 1) * data[i, 1, _zeros[1] - 1]
                )
                data[i, 2, _zeros[1:]] = (
                    np.ones(len(_zeros) - 1) * data[i, 2, _zeros[1] - 1]
                )
                data[i, 3, _zeros[1:]] = (
                    np.ones(len(_zeros) - 1) * data[i, 3, _zeros[1] - 1]
                )
                data[i, 4, _zeros[1:]] = (
                    np.ones(len(_zeros) - 1) * data[i, 4, _zeros[1] - 1]
                )
                data[i, 5, _zeros[1:]] = (
                    np.ones(len(_zeros) - 1) * data[i, 5, _zeros[1] - 1]
                )
                data[i, 6, _zeros[1:]] = (
                    np.ones(len(_zeros) - 1) * data[i, 6, _zeros[1] - 1]
                )
                data[i, 7, _zeros[1:]] = (
                    np.ones(len(_zeros) - 1) * data[i, 7, _zeros[1] - 1]
                )
                data[i, 8, _zeros[1:]] = (
                    np.ones(len(_zeros) - 1) * data[i, 8, _zeros[1] - 1]
                )

        if ("binaries_" + self.metstr + ".h5") not in os.listdir(
            self.output_dir
        ) or self.overwrite:
            print("Saving processed data to", self.output_dir)

            ds = xr.Dataset(
                data_vars=dict(
                    time=(["model", "step"], data[:, 0, :]),
                    mass1=(["model", "step"], data[:, 1, :]),
                    log_Lbol1=(["model", "step"], data[:, 2, :]),
                    log_Teff1=(["model", "step"], data[:, 3, :]),
                    log_R1=(["model", "step"], data[:, 4, :]),
                    mass2=(["model", "step"], data[:, 5, :]),
                    log_Lbol2=(["model", "step"], data[:, 6, :]),
                    log_Teff2=(["model", "step"], data[:, 7, :]),
                    log_R2=(["model", "step"], data[:, 8, :]),
                ),
                coords=dict(
                    model=(["model"], model_orbits),
                    step=("step", np.arange(first_empty)),
                ),
                attrs=dict(
                    description="BPASS evolution data for binary stars at Z="
                    + str(self.met)
                ),
            )

            ds.to_netcdf(self.output_dir + "/binaries_" + self.metstr + ".h5")

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
            convert_single_data     Processes single stellar track data into a SingleStarTable
            convert_binary_data     Processes binary stellar track data into a BinaryStarTable
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
        SingleStarTable.
        """

        # Create directory if it does not already exists
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

        model_directory = self.input_dir + self.met + "/single/"

        # Create the zero array; assume 1000 models and clean up after
        num_files = int(1000)
        # Times for time array
        num_times = int(1e5)  # Assume there are at most 1e4 outputs
        # Save properties as a function of time: (1) time, (2) mass,
        # (3) bolometric luminosity, (4) surface temperature, (5) radius
        data = np.zeros((num_files, 5, num_times))

        i = 0
        for model in os.listdir(model_directory):

            if model.endswith("_0_compressed.pkl.gz"):

                df = pd.read_pickle(model_directory + "/" + model, compression="gzip")

                # We want to extract the following values
                _num = np.arange(len(df.star_age.values))

                data[i, 0, _num] = df.star_age.values.astype("float")  # time in yr
                data[i, 1, _num] = df.star_mass.values.astype("float")  # mass in Msun
                data[i, 2, _num] = df.log_L.values.astype("float")  # log Lbol in Lsun
                data[i, 3, _num] = df.log_Teff.values.astype("float")  # log Teff in K
                data[i, 4, _num] = df.log_R.values.astype("float")  # log R in Rsun

                i += 1

        # Remove superfluous models
        first_empty = 0
        i = 0
        while i < num_files:
            if len(np.nonzero(data[i, 1, :])[0]) > 0:
                first_empty += 1
                i += 1
            else:
                first_empty += 1
                i = num_times

        data = data[:first_empty, :, :]

        # Remove superfluous model numbers
        first_empty = 0
        j = 0
        while j < num_times:
            if len(np.nonzero(data[:, 1, j])[0]) > 0:
                first_empty += 1
                j += 1
            else:
                first_empty += 1
                j = num_times

        data = data[:, :, :first_empty]

        # If fewer models, fill with sensible values
        for i in range(len(data[:, 0, 0])):
            _zeros = np.where(data[i, 0, :] == 0)[0]
            if len(_zeros) > 1:  # if not longest
                data[i, 0, _zeros[1:]] = (
                    np.ones(len(_zeros) - 1) * data[i, 0, _zeros[1] - 1]
                )
                data[i, 1, _zeros[1:]] = (
                    np.ones(len(_zeros) - 1) * data[i, 1, _zeros[1] - 1]
                )
                data[i, 2, _zeros[1:]] = (
                    np.ones(len(_zeros) - 1) * data[i, 2, _zeros[1] - 1]
                )
                data[i, 3, _zeros[1:]] = (
                    np.ones(len(_zeros) - 1) * data[i, 3, _zeros[1] - 1]
                )
                data[i, 4, _zeros[1:]] = (
                    np.ones(len(_zeros) - 1) * data[i, 4, _zeros[1] - 1]
                )

        # Sort by mass
        _sort = np.argsort(data[:, 1, 0])
        data = data[_sort, :, :]

        if ("singles_" + self.met + ".h5") not in os.listdir(
            self.output_dir
        ) or self.overwrite:
            print("Saving processed data to", self.output_dir)

            masses = data[:, 1, 0]

            ds = xr.Dataset(
                data_vars=dict(
                    time=(["model", "step"], data[:, 0, :]),
                    mass=(["model", "step"], data[:, 1, :]),
                    log_Lbol=(["model", "step"], data[:, 2, :]),
                    log_Teff=(["model", "step"], data[:, 3, :]),
                    log_R=(["model", "step"], data[:, 4, :]),
                ),
                coords=dict(
                    model=("model", masses), step=("step", np.arange(first_empty))
                ),
                attrs=dict(
                    description="MPA/Bonn evolution data for single stars at Z="
                    + str(self.met)
                ),
            )
            ds.to_netcdf(self.output_dir + "/singles_" + self.met + ".h5")

        else:
            print("Cannot save model. Try setting overwrite=True...")

        return

    def convert_binary_data(self):
        pass
