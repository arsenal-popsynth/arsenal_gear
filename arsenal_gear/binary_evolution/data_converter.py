"""
data_converter.py
================================

This file defines the interface to various binary evolution models
through downloading, reorganzing and interpreting their outputs.
"""

import os
import tarfile
from zipfile import ZipFile
from pathlib import Path
import requests
from tqdm import tqdm

import astropy.units as u
import numpy as np
import xarray as xr

from astropy.units import Quantity
from arsenal_gear.population import BinaryPopulation, StarPopulation

def convert_singles(BPASS_directory, output_directory='./arsenal_BPASS', metals='z014', overwrite=False):

    # Create directory if it does not already exists
    Path(output_directory).mkdir(parents=True, exist_ok=True)

    model_directory = BPASS_directory + '/singles/'

    for metal_directory in os.listdir(model_directory):
        
        if not metal_directory.startswith("z"):
            continue

        # Create the zero array
        # Number of files depends on the choice of mass limits for the IMF
        num_files = len(os.listdir(model_directory + '/' + metal_directory)) # -1 for directory itself
        # Times for time array
        num_times = 502 # from log(age/yr) = 6 to 9, with 10 times more models than the public models and a 0
        times = np.concatenate((np.zeros(1), np.logspace(4, 9, 501)))
        # Save 4 properties: mass, luminosity, temperature, radius
        data = np.zeros((num_files, 3, num_times))

        i = 0
        for model in os.listdir(model_directory + '/' + metal_directory):

            if model.startswith("sneplot"):

                _data = np.genfromtxt(model_directory + '/' + metal_directory 
                                  + '/' + model)

                # Replace NaNs by 0s --> What about files without a companion?
                _data = np.nan_to_num(_data)
                # We want to extract the following values
                _times = _data[:, 1]
                for t in range(len(_times)):

                    if t == 0:
                        _mask = np.where(_times[t] >= times)[0]
                    elif (t == (len(_times) - 1)) and (_times[-1] < times[-1]):
                        _mask = np.arange(len(times))[np.where(_times[t] >= times)[0][-1]:]
                    else:
                        _mask = np.where((_times[t] >= times) & (_times[t-1] < times))[0]

                    for m in _mask:

                        data[i, 0, m] = _data[t, 5]     # mass in MSun
                        data[i, 1, m] = _data[t, 4] # Lsun
                        if t == (len(_times) - 1):
                            data[i, 1, m] *= 0 # must set to 0 after SN
                        data[i, 2, m] = _data[t, 3] # K

                i += 1

        # Sort by mass
        _sort = np.argsort(data[:, 0, 0])
        data = data[_sort, :, :]


        if ('singles_' + metals + '.h5') not in os.listdir(output_directory) or overwrite:
            print("Saving processed data to", output_directory)

            # Times
            times = np.round(times, 2).astype('str')
            for t in range(len(times)):
                times[t] = times[t].ljust(4, '0')

            # Masses as strings
            masses = data[:, 0, 0].astype('str')
            for m in range(len(masses)):
                masses[m] = masses[m].ljust(4, '0')

            ds = xr.DataArray(data, coords=[("Model", masses), ("Property", ["Mass (MSun)", "log L_bol (LSun)", "log T_eff (K)"]), ("Time (log t/yr)", times)])
            ds.to_netcdf(output_directory + '/singles_' + metals + '.h5')

        else:
            print("Cannot save model. Try setting overwrite=True...")

    return

def convert_binaries(BPASS_directory, output_directory='./arsenal_BPASS', metals='z014', overwrite=False):

    # Create directory if it does not already exists
    Path(output_directory).mkdir(parents=True, exist_ok=True)

    model_directory = BPASS_directory + '/binaries/'

    # Values for disrupted or modified systems, when looking for new systems
    M_single = np.concatenate((np.arange(0.1, 10., 0.1) , 
                               np.arange(10, 100, 1),
                               np.arange(100, 325, 25)))

    for metal_directory in os.listdir(model_directory):
        
        if not metal_directory.startswith("z"):
            continue

        # Create the zero array
        # Number of files depends on the choice of IMF
        num_files = len(os.listdir(model_directory + '/' + metal_directory)) # -1 for directory itself
        # Times for time array
        num_times = 502 # from log(age/yr) = 4 to 9, with 10 times more models than the public models and a 0
        times = np.concatenate((np.zeros(1), np.logspace(4, 9, 501)))
        # Save 7 properties: masses, luminosities, temperatures, period
        data = np.zeros((num_files, 7, num_times))

        M2_for_models = np.zeros(num_files)

        i = 0
        for model in os.listdir(model_directory + '/' + metal_directory):

            if model.startswith("sneplot"):

                merger = False
                supernova = False
                rejuvenated = False
                accreted_mass = 0
                accretion_time = None
                effective_time = None
                rejuvenated_file = None

                _data = np.genfromtxt(model_directory + metal_directory 
                                  + '/' + model)

                # Replace NaNs by 0s --> What about files without a companion?
                _data = np.nan_to_num(_data)

                # Effective companion mass
                M2_eff = _data[0, 37]
                M2_for_models[i] = M2_eff
                
                # We want to extract the following values
                _times = _data[:, 1]
                for t in range(len(_times)):

                    if t == 0:
                        _mask = np.where(_times[t] >= times)[0]
                    elif (t == (len(_times) - 1)) and (_times[-1] < times[-1]):
                        _mask = np.arange(len(times))[np.where(_times[t] >= times)[0][-1]:]
                        # If initial mass above 8 Msun, assume SN - to adjust? 
                        if _data[0, 5] > 8:
                            supernova = True
                    else:
                        _mask = np.where((_times[t] >= times) & (_times[t-1] < times))[0]
                        if (_data[t, 5] > _data[t-1, 5]) and (_data[t, 37] == _data[t-1, 37]):# and (merger == False):
                            # Assume merger if M1 increased but M2 remained fixed
                            # Note that this can take place over several timesteps
                            merger = True
                            accreted_mass += _data[t, 5] - _data[t-1, 5]
                        elif (_data[t, 37] > _data[t-1, 37]) and (merger == False):
                            # Assume accretion if M2 increased
                            M2_eff = _data[t, 37]
                            if (_data[t, 37]/_data[0, 37] >= 1.05) and (_data[0, 37] >= 2):
                                rejuvenated = True
                            if accretion_time == None:
                                accretion_time = _times[t]

                    if supernova and not merger:

                        # Set mass to use for post-SN star
                        M2_closest = M_single[np.argmin(np.abs(M_single - M2_eff))]
                        if np.round(M2_closest, 1) == np.round(M2_closest, 0):
                            M2_closest = str(int(M2_closest))
                        elif M2_closest < 10:
                            M2_closest = str(np.round(M2_closest, 1))
                        else:
                            M2_closest = str(int(M2_closest))

                        # Set current time for the rejuvenated star
                        if rejuvenated:
                            effective_time = _times[t] - accretion_time
                        else:
                            effective_time = _times[t]

                        rejuvenated_file = BPASS_directory + '/singles/' + metal_directory + '/sneplot-' + metal_directory + '-' + M2_closest


                    else:
                        
                        for m in _mask:

                            data[i, 0, m] = _data[t, 5]      # mass in MSun
                            data[i, 1, m] = _data[t, 4]  # Lsun
                            if t == (len(_times) - 1):
                                data[i, 1, m] *= 0 # must set to 0 after SN
                            data[i, 2, m] = _data[t, 3]  # K 
                            # Companion
                            if not merger:
                                data[i, 3, m] = _data[t, 37]     # companion mass in MSun
                                data[i, 4, m] = _data[t, 48] # Lsun, companion
                                data[i, 5, m] = _data[t, 47] # K, companion
                            if merger: 
                                if (_data[t, 5] - _data[t-1, 5]) > 0: # If still accreting
                                    data[i, 3, m] = _data[t, 37] - accreted_mass
                                    data[i, 4, m] = _data[t, 48] # Lsun, companion
                                    data[i, 5, m] = _data[t, 47] # K, companion
                            data[i, 6, m] = _data[t, 34]  # yr


                    if rejuvenated_file:

                        _data_M2 = np.genfromtxt(rejuvenated_file)

                        # Replace NaNs by 0s --> What about files without a companion?
                        _data_M2 = np.nan_to_num(_data_M2)

                        _times_M2 = _data_M2[:, 1][_data_M2[:, 1] > effective_time]
                        
                        for t in range(len(_times_M2)):

                            if (t == (len(_times_M2) - 1)) and (_times_M2[-1] < times[-1]):
                                _mask = np.arange(len(times))[np.where(_times_M2[t] >= times)[0][-1]:]
                            else:
                                _mask = np.where((_times_M2[t] >= times) & (_times_M2[t-1] < times))[0]


                            for m in _mask:

                                # Keep primary mass
                                data[i, 0, m] = data[i, 0, m-1] # mass in MSun
                                # Companion properties
                                data[i, 3, m] = _data_M2[t, 5]  # mass in MSun
                                data[i, 4, m] = _data_M2[t, 4]  # Lsun
                                if t == (len(_times_M2) - 1):
                                    data[i, 4, m] *= 0 # must set to 0 after SN
                                data[i, 5, m] = _data_M2[t, 3]  # K 


                i += 1
        

        if ('binaries_' + metals + '.h5') not in os.listdir(output_directory) or overwrite:
            print("Saving processed data to", output_directory)

            # Times
            times = np.round(times, 2).astype('str')
            for t in range(len(times)):
                times[t] = times[t].ljust(4, '0')

            # Masses as strings
            masses = data[:, 0, 0].astype('str')
            for m in range(len(masses)):
                masses[m] = masses[m].ljust(4, '0')

            mass_ratios = np.round(M2_for_models/data[:, 0, 0], 1).astype('str')

            periods = np.round(np.log10(data[:, 6, 0]*365.25), 1).astype('str')

            models = np.zeros(len(masses)).astype('str')
            for i in range(len(masses)):
                models[i] = masses[i] + '-' + mass_ratios[i] + '-' + periods[i]

            _sort = np.argsort(models)
            data_to_save = data[_sort, :, :]

            ds = xr.DataArray(data_to_save, coords=[("Model", models[_sort]), ("Property", ["Mass 1 (MSun)", "log L_bol 1 (LSun)", "log T_eff 1 (K)",
                                            "Mass 2 (MSun)", "log L_bol 2 (LSun)", "log T_eff 2 (K)", "P (yr)"]), ("Time (log t/yr)", times)])
            ds.to_netcdf(output_directory + '/binaries_' + metals + '.h5')

        else:
            print("Cannot save model. Try setting overwrite=True...")

    return


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

        