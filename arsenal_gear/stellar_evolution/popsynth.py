"""
output from population synthesis
================================

This submodule defines the interface to various binary population
synthesis codes through interpreting and processing their outputs
"""

import os
import tarfile
from zipfile import ZipFile

import astropy.units as u
import numpy as np
import requests
from astropy.units import Quantity
from tqdm import tqdm

from arsenal_gear.population import BinaryPopulation, StarPopulation
from arsenal_gear.utils.convert_BPASS_output import *

class BPASS_stellar_models:
    """
    Reads in BPASS stellar model files for use with a discrete stellar population.
    """

    # basic options for BPASS stellar models
    # bpass_url = "https://www.dropbox.com/scl/fo/mpuas1xh5owmdadu0vpev/h?dl=0" + \
    #            "&e=1&rlkey=7vlk7ra6kvoztzmae8wr34kmz"
    # Change dl=0 to dl=1 to force download
    bpass_url = (
        "https://www.dropbox.com/scl/fo/mpuas1xh5owmdadu0vpev/h?dl=1"
        + "&e=1&rlkey=7vlk7ra6kvoztzmae8wr34kmz"
    )

    def __init__(
        self,
        singles: StarPopulation,
        binaries: BinaryPopulation,
        metal: str,
        time: Quantity["time"],
        bpass_dir: str,
        force_download: bool = False,
    ) -> None:
        """
        Args:
            stars:     the StarPopulation instance for which we want information
            binaries:  the BinaryPopulation instance for which we want information
            metal:     the metallicity of the stellar population
            TO DO -CCC, 26/06/2025: Read metallicity from binaries
            time:      current simulation time
            bpass_dir: the directory for the BPASS models

        """

        self.download = force_download

        self.metals = [
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
        if metal not in self.metals:
            raise ValueError(f"Metallicity {metal} does not exist.")

        self.metal: str = metal
        if bpass_dir[-1] == "/":
            self.dir: str = bpass_dir
        else:
            self.dir: str = bpass_dir + "/"

        self.time: Quantity["time"] = time

        self.s_mass: list = singles["mass"].to(u.Msun).value
        self.b_mass: list = binaries.primary["mass"].to(u.Msun).value
        self.mratio: list = binaries.secondary["mass"] / binaries.primary["mass"]
        self.logp: list = np.log10(binaries["period"].to(u.d).value)
        super().__init__()

    def downloader(
        self, url=bpass_url, message="Downloading stellar models..."
    ) -> None:
        """
        Method for downloading BPASS stellar models from the web.

        Args:
            url (str): The URL of the BPASS repository.
            message (str): Optional message to display before downloading.

        Raises:
            Exception: If the download fails.

        """
        if message is not None:
            print(message)

        download_url = url
        if self.dir[-1] == "/":
            fname = self.dir + "bpass_v2.2.zip"
        else:
            fname = self.dir + "/bpass_v2.2.zip"

        try:
            response = requests.get(download_url, stream=True, timeout=10)
        except requests.exceptions.Timeout as e:
            raise requests.exceptions.Timeout(
                f"Request timed out. Check your internet connection. {e}"
            )
        except requests.exceptions.RequestException as e:
            raise requests.exceptions.RequestException(f"Request failed: {e}")

        # Get file size
        total_size = int(response.headers.get("content-length", 0))
        # create a progress bar
        tqdm_args = {'desc': 'Downloading', 'total': total_size, 'unit': 'B',
                        'unit_scale': True, 'unit_divisor': 1024}
        # write the file
        with open(fname, "wb") as f, tqdm(**tqdm_args) as prog_bar:
            for chunk in response.iter_content(chunk_size=1024):
                f.write(chunk)
                prog_bar.update(len(chunk))

    def unzip(
        self, zip_name: str, target_file=None, delete_zip=False, inspect=False
    ) -> None:
        """
        Un-compress the zip file downloaded from the BPASS dropbox.

        Args:
            zip_name    (str): Name of expected zip file
            target_file (str): Name of expected tar file in zip file
                               If none, extract all tar files
            delete_zip (bool): Delete the zip file after extracting the tar file
            inspect    (bool): Print the names of the tar files in the zip file
        """
        if self.dir[-1] == "/":
            fname = self.dir + zip_name
        else:
            fname = self.dir + "/" + zip_name

        with ZipFile(fname, "r") as zip_archive:
            file_names = zip_archive.namelist()
            if inspect:
                for name in file_names:
                    print(name)
            else:
                if target_file in file_names:

                    zip_archive.extract(target_file, path=self.dir)
                elif target_file is None:
                    zip_archive.extractall(path=self.dir)

        if delete_zip and fname.exists():
            print("Deleting", fname)
            fname.unlink()

    def untar(self, tar_name: str, delete_tar=False) -> None:
        """
        Un-compress the tar file extracted from the BPASS zip file.

        Args:
            tar_name    (str): Name of expected tar file
            delete_tar (bool): Delete the tar file after extracting
        """
        if self.dir[-1] == "/":
            fname = self.dir + tar_name
        else:
            fname = self.dir + "/" + tar_name

        if not tarfile.is_tarfile(fname):
            raise OSError(
                f"{fname} is not a valid tar.gz file. "
                "Try again with `force_download=True`"
            )
        with tarfile.open(fname, "r:gz") as tar:
            print(f"Extracting {fname}...")
            tar.extractall(path=self.dir)

        if delete_tar and fname.exists():
            print("Deleting", fname)
            fname.unlink()

    def data_from_model(self, data: tuple, time_now: Quantity["time"]) -> tuple:
        """
        Extract the stellar properties necessary for feedback
        from a BPASS stellar model file.

        Args:
            data                (tuple): Matrix loaded from the stellar model file
            time_now (Quantity["time"]): Current simulation time

        Returns:
            Feedback properties
        """
        # Replace NaNs by 0s --> What about files without a companion?
        data = np.nan_to_num(data)
        # We want to extract the following values
        time = data[:, 1] * u.yr
        ind_now = np.where((time_now - time) > 0 * u.yr)[-1][-1]

        # Model time
        t_max = time[-1]

        # Need to use column 36 for 50-0.9-1
        M1 = (
            data[ind_now, 36] * u.Msun
        )  # 5 recommended in manual but doesn't conserve mass, 36 works
        L1 = 10 ** data[ind_now, 4] * u.Lsun
        T1 = 10 ** data[ind_now, 3] * u.K
        R1 = 10 ** data[ind_now, 2] * u.Rsun

        M2 = data[ind_now, 37] * u.Msun
        L2 = 10 ** data[ind_now, 48] * u.Lsun
        T2 = 10 ** data[ind_now, 47] * u.K
        R2 = 10 ** data[ind_now, 46] * u.Rsun

        M2_init = data[0, 37] * u.Msun  # Initial companion mass
        P = data[ind_now, 34] * u.yr  # Final orbital period

        dM = (
            (
                data[ind_now, 39]
                + data[ind_now, 40]
                - data[ind_now, 41]
                - data[ind_now, 42]
                + data[ind_now, 43]
                + data[ind_now, 44]
            )
            * u.Msun
            / (1.989 * u.s)
        )

        return M1, M2, M2_init, P, R1, R2, T1, T2, L1, L2, dM, t_max

    def data_from_single_model(self, data: tuple, time_now: Quantity["time"]) -> tuple:
        """
        Extract the stellar properties necessary for feedback
        from a BPASS stellar model file, for a single star.

        Args:
            data                (tuple): Matrix loaded from the stellar model file
            time_now (Quantity["time"]): Current simulation time

        Returns:
            Feedback properties
        """
        # Replace NaNs by 0s --> What about files without a companion?
        data = np.nan_to_num(data)
        # We want to extract the following values
        time = data[:, 1] * u.yr
        ind_now = np.where((time_now - time) > 0 * u.yr)[-1][-1]

        # Model time
        t_max = time[-1]

        # Need to use column 36 for 50-0.9-1
        M1 = (
            data[ind_now, 5] * u.Msun
        )  # 5 recommended in manual but doesn't conserve mass, 36 works
        L1 = 10 ** data[ind_now, 4] * u.Lsun
        T1 = 10 ** data[ind_now, 3] * u.K
        R1 = 10 ** data[ind_now, 2] * u.Rsun

        dM = data[ind_now, 39] * u.Msun / (1.989 * u.s)

        return M1, R1, T1, L1, dM, t_max

    def data_from_companion_model(
        self, data: tuple, time_now: Quantity["time"]
    ) -> tuple:
        """
        Extract the stellar properties necessary for feedback
        from a BPASS stellar model file.

        Args:
            data                (tuple): Matrix loaded from the stellar model file
            time_now (Quantity["time"]): Current simulation time

        Returns:
            Feedback properties
        """
        # Replace NaNs by 0s --> What about files without a companion?
        data = np.nan_to_num(data)
        # We want to extract the following values
        time = data[:, 1] * u.yr
        ind_now = np.where((time_now - time) > 0 * u.yr)[-1][-1]

        # Model time
        t_max = time[-1]

        # Directory
        c_dir = self.dir + "NEWBINMODS/NEWSECMODS/" + self.metal + "_2/"

        # Need to use column 36 for 50-0.9-1
        MR = (
            data[ind_now, 29] * u.Msun
        )  # 5 recommended in manual but doesn't conserve mass, 36 works

        M2 = data[ind_now, 37] * u.Msun

        M2_init = data[0, 37] * u.Msun  # Initial companion mass
        P = data[ind_now, 34] * u.yr  # Final orbital period

        M2_vals = np.concatenate(
            (
                np.arange(0.1, 0.7, 0.1),
                np.arange(0.8, 2.2, 0.1),
                np.array([2.3, 2.5, 2.7, 3, 3.2, 3.5, 3.7]),
                np.arange(4, 10, 0.5),
                np.arange(10, 26, 1),
                np.array([30, 35, 40, 50, 60, 70, 80, 100, 120, 150, 200, 300]),
            )
        )
        MR_vals = np.array(
            [
                0.1,
                0.158489,
                0.199526,
                0.251189,
                0.316228,
                0.398107,
                0.501187,
                0.630957,
                0.794328,
                1.29593,
                1.4,
                1.99526,
                2.51189,
                3.16228,
                3.98107,
                5.01187,
                6.30957,
                7.94328,
                10.0000,
                12.5893,
            ]
        )

        P_vals = np.arange(-0.4, 4.2, 0.2)

        M2_file = M2_vals[np.argmin(np.abs(M2_vals - M2_init.to_value(u.Msun)))]
        if M2_file == int(M2_file):
            M2_str = str(int(M2_file))
        else:
            M2_str = str(M2_file)
        MR_file = MR_vals[np.argmin(np.abs(MR_vals - MR.to_value(u.Msun)))]
        if MR_file == 10:
            MR_str = "10.0000"  # set manually
        else:
            MR_str = str(MR_file)
        P_file = P_vals[np.argmin(np.abs(P_vals - np.log10(P.to_value(u.day))))]
        P_str = str(P_file)
        if (P_str[-1] != 0) and (len(P_str) > 3):
            P_str = P_str[:3]

        c_fname = "sneplot_2-" + self.metal + "-" + M2_str + "-" + MR_str + "-" + P_str

        if c_fname in os.listdir(c_dir):
            c_data = np.genfromtxt(c_dir + c_fname)

        else:
            num0 = 0
            print("Adding 0's...")
            while num0 < 10:
                c_fname = c_fname + "0"
                if (c_fname) in os.listdir(c_dir):
                    c_data = np.genfromtxt(c_dir + c_fname)
                    num0 += 10  # force exit the loop if ok
                else:
                    num0 += 1

            if (c_fname) not in os.listdir(c_dir):
                print("Incorrect file!", c_fname)
                c_fname = c_fname[:-10]
                try_num = [-1, 1, -2, 2, -3, 3, -4, 4, -5, 5]
                for n in try_num:
                    new_num = str(float(c_fname[-3:]) + 0.2 * n)
                    if len(new_num) > 3:
                        new_num = new_num[:3]
                    n_fname = c_fname[:-3] + new_num

                    if n_fname in os.listdir(c_dir):
                        c_data = np.genfromtxt(c_dir + n_fname)

                    else:
                        num0 = 0
                        print("Adding 0's...")
                        while num0 < 10:
                            n_fname = n_fname + "0"
                            if (n_fname) in os.listdir(c_dir):
                                c_data = np.genfromtxt(c_dir + n_fname)
                                num0 += 10  # force exit the loop if ok
                            else:
                                num0 += 1

        # Replace NaNs by 0s --> What about files without a companion?
        c_data = np.nan_to_num(c_data)
        # We want to extract the following values
        time = c_data[:, 1] * u.yr
        ind_now = np.where((time_now - time) > 0 * u.yr)[-1][-1]

        # Model time
        t_max = time[-1]

        # Need to use column 36 for 50-0.9-1
        M2 = (
            c_data[ind_now, 5] * u.Msun
        )  # 5 recommended in manual but doesn't conserve mass, 36 works
        L2 = 10 ** c_data[ind_now, 4] * u.Lsun
        T2 = 10 ** c_data[ind_now, 3] * u.K
        R2 = 10 ** c_data[ind_now, 2] * u.Rsun

        dM = c_data[ind_now, 39] * u.Msun / (1.989 * u.s)

        return M2, R2, T2, L2, dM, t_max

    def read_bpass_data(self) -> tuple:
        """
        Read BPASS data from a stellar model file.
        """
        # General file name, to append
        b_fname = (
            self.dir + "NEWBINMODS/NEWBINMODS/" + self.metal + "/sneplot-" + self.metal
        )
        s_fname = self.dir + "NEWSINMODS/" + self.metal + "/sneplot-" + self.metal
        # Get frequency for each unique binary and single
        b_props = np.vstack((self.b_mass, self.mratio, self.logp))
        b_vals, _, b_counts = np.unique(
            b_props, return_index=True, return_counts=True, axis=1
        )
        s_vals, _, s_counts = np.unique(
            self.s_mass, return_index=True, return_counts=True
        )
        len_b = len(b_counts)
        len_s = len(s_counts)
        feedback = np.empty((10, len_b + len_s))
        feedback[0, :] = np.concatenate((b_counts, s_counts))

        # Get the properties for each unique binary
        for i in range(len_b + len_s):

            if i in range(len_b):
                _mass, _mrat, _logp = np.round(b_vals[:, i], decimals=1)
                # Matching file names from BPASS
                if str(_mass)[-2:] == ".0":
                    _mass = int(_mass)
                if str(_logp)[-2:] == ".0":
                    _logp = int(_logp)
                fname = b_fname + "-" + str(_mass) + "-" + str(_mrat) + "-" + str(_logp)

                data = self.data_from_model(np.genfromtxt(fname), self.time)
                feedback[1, i] = data[0].to_value(u.Msun)
                feedback[2, i] = data[1].to_value(u.Msun)
                feedback[3, i] = data[4].to_value(u.Rsun)
                feedback[4, i] = data[5].to_value(u.Rsun)
                feedback[5, i] = data[6].to_value(u.K)
                feedback[6, i] = data[7].to_value(u.K)
                feedback[7, i] = data[8].to_value(u.Lsun)
                feedback[8, i] = data[9].to_value(u.Lsun)
                feedback[9, i] = data[10].to_value(u.Msun / u.yr)

                # Set primary properties to 0 if remnant
                if (data[-1] > 0 * u.yr) and (data[-1] < self.time):

                    print("Getting data from companion... assuming no kick!")
                    feedback[7, i] = 0

                    data_rem = self.data_from_companion_model(
                        np.genfromtxt(fname), self.time
                    )
                    feedback[2, i] = data_rem[0].to_value(u.Msun)
                    feedback[4, i] = data_rem[1].to_value(u.Rsun)
                    feedback[6, i] = data_rem[2].to_value(u.K)
                    feedback[8, i] = data_rem[3].to_value(u.Lsun)
                    feedback[9, i] = data_rem[4].to_value(u.Msun / u.yr)

            else:
                _mass = np.round(s_vals[i - len_b], decimals=1)
                # Matching file names from BPASS
                if str(_mass)[-2:] == ".0":
                    _mass = int(_mass)
                if _mass == 120:
                    _mass = 125  # Correct manually; not same mass as the binary models
                fname = s_fname + "-" + str(_mass)

                data = self.data_from_single_model(np.genfromtxt(fname), self.time)
                feedback[1, i] = data[0].to_value(u.Msun)
                feedback[3, i] = data[1].to_value(u.Rsun)
                feedback[5, i] = data[2].to_value(u.K)
                feedback[7, i] = data[3].to_value(u.Lsun)
                feedback[9, i] = data[4].to_value(u.Msun / u.yr)
                feedback[2, i] = 0  # Set companion properties to 0
                feedback[4, i] = 0
                feedback[6, i] = 0
                feedback[8, i] = 0

                # Also set primary properties to 0 if remnant
                if (data[-1] > 0 * u.yr) and (data[-1] < self.time):
                    # Setting properties to 0 for now
                    feedback[7, i] = 0
                    feedback[9, i] = 0

        return feedback

    def get_stellar_models(
        self,
        tar_name="bpass-v2.2-newmodels.tar.gz",
        zip_name="bpass_v2.2.zip",
        url=bpass_url,
    ) -> tuple:
        """
        Look for the stellar models for the stellar population
        and download/unzip/untar them if they are not available.
        """
        if os.path.isdir(self.dir + "NEWBINMODS"):
            print("Reading BPASS data...")
            _tmp = self.read_bpass_data()
        else:
            print("BPASS data not available at", self.dir + "NEWBINMODS")
            print("Looking for a tar file...")
            if os.path.isfile(self.dir + tar_name):
                self.untar(tar_name, delete_tar=False)
                _tmp = self.read_bpass_data()

            else:
                print("tar file not available at", self.dir)
                print("Looking for a zip file...")
                if os.path.isfile(self.dir + zip_name):
                    self.unzip(
                        zip_name, target_file=tar_name, delete_zip=False, inspect=False
                    )
                    self.untar(tar_name, delete_tar=False)
                    _tmp = self.read_bpass_data()

                else:
                    print("zip file not available at", self.dir)
                    if self.download:
                        self.downloader(url)
                        self.unzip(
                            zip_name,
                            target_file=tar_name,
                            delete_zip=False,
                            inspect=False,
                        )
                        self.untar(tar_name, delete_tar=False)
                        _tmp = self.read_bpass_data()

                    else:
                        print("Set force_download = True to download the files")
                        _tmp = np.empty((10, 1))

        return (
            _tmp[0],
            _tmp[1] * u.Msun,
            _tmp[2] * u.Msun,
            _tmp[3] * u.Rsun,
            _tmp[4] * u.Rsun,
            _tmp[5] * u.K,
            _tmp[6] * u.K,
            _tmp[7] * u.Lsun,
            _tmp[8] * u.Lsun,
            _tmp[9] * u.Msun / u.yr,
        )
