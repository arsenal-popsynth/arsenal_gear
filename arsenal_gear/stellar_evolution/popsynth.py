"""
output from population synthesis
==========

This submodule defines the interface to various binary population
synthesis codes through interpreting and processing their outputs
"""

import os
import tarfile
from zipfile import ZipFile

import numpy as np
import requests
from tqdm import tqdm

import astropy.units as u
from astropy.units import Quantity

from arsenal_gear.population import StarPopulation, BinaryPopulation
    
class BPASS_stellar_models():
    """
    Reads in BPASS stellar model files for use with a discrete stellar population.
    """
    # basic options for BPASS stellar models
    #bpass_url = "https://www.dropbox.com/scl/fo/mpuas1xh5owmdadu0vpev/h?dl=0&e=1&rlkey=7vlk7ra6kvoztzmae8wr34kmz"
    # Change dl=0 to dl=1 to force download
    bpass_url = "https://www.dropbox.com/scl/fo/mpuas1xh5owmdadu0vpev/h?dl=1&e=1&rlkey=7vlk7ra6kvoztzmae8wr34kmz"
    
    def __init__(self, binaries: BinaryPopulation, metal: str,
                 time: Quantity["time"],
                 bpass_dir: str, force_download: bool = False) -> None:
        """
        Args:
            binaries:  the BinaryPopulation instance for which we want information
            metal:     the metallicity of the stellar population
            TO DO -CCC, 26/06/2025: Read metallicity from binaries
            time:      current simulation time
            #prev_time: previous simulation time, necessary for binary mass loss
            bpass_dir: the directory for the BPASS models

        """
        
        self.download = force_download
            
        self.metals = ["zem5", "zem4", "z001", "z002", "z003", "z004", "z006",
                       "z008", "z010", "z014",  "z020",  "z030",  "z040"]
        if metal not in self.metals:
            raise ValueError(f"Metallicity {metal} does not exist.")
            
        self.metal:    str = metal
        if bpass_dir[-1] == '/':
            self.dir:    str = bpass_dir
        else:
            self.dir:    str = bpass_dir + '/'
            
        self.time:    Quantity["time"] = time
            
        self.mass:    list = binaries.primary["mass"].to(u.Msun).value
        self.mrat:    list = binaries.secondary["mass"] / binaries.primary["mass"]
        self.logp:    list = np.log10(binaries["period"].to(u.d).value)
        super().__init__()
                

    def downloader(self, url = bpass_url, 
                   message = "Downloading stellar models...") -> None:
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
        if self.dir[-1] == '/':
            fname = self.dir + 'bpass_v2.2.tar.gz'
        else:
            fname = self.dir + '/bpass_v2.2.tar.gz'
            
        try:
            response = requests.get(download_url, stream=True, timeout=10)
        except requests.exceptions.Timeout:
            raise Exception('Request timed out. Check your internet connection.')
        except requests.exceptions.RequestException as e:
            raise Exception(f'Request failed: {e}')

        # Get file size
        total_size = int(response.headers.get('content-length', 0))
        # create a progress bar
        tqdm_args = dict(desc='Downloading', total=total_size, unit='B',
                        unit_scale=True, unit_divisor=1024)
        # write the file
        with open(fname, 'wb') as f, tqdm(**tqdm_args) as prog_bar:
            for chunk in response.iter_content(chunk_size=1024):
                f.write(chunk)
                prog_bar.update(len(chunk))
                
    def unzip(self, zip_name: str, target_file=None, 
              delete_zip=False, inspect=False) -> None:
        """
        Un-compress the zip file downloaded from the BPASS dropbox.
        
        Args:
            zip_name    (str): Name of expected zip file
            target_file (str): Name of expected tar file in zip file
                               If none, extract all tar files
            delete_zip (bool): Delete the zip file after extracting the tar file
            inspect    (bool): Print the names of the tar files in the zip file
        """
        if self.dir[-1] == '/':
            fname = self.dir + zip_name
        else:
            fname = self.dir + '/' + zip_name
                        
        with ZipFile(fname, 'r') as zip_archive:
            file_names = zip_archive.namelist()
            if inspect:
                for name in file_names:
                    print(name)
            else:
                if target_file in file_names:
                    
                    zip_archive.extract(target_file, path = self.dir)
                elif target_file is None:
                    zip_archive.extractall(path = self.dir)

        if delete_zip and fname.exists():
            print("Deleting", fname)
            fname.unlink()
                
    
    def untar(self, tar_name: str, delete_tar = False) -> None:
        """
        Un-compress the tar file extracted from the BPASS zip file.
        
        Args:
            tar_name    (str): Name of expected tar file
            delete_tar (bool): Delete the tar file after extracting
        """
        if self.dir[-1] == '/':
            fname = self.dir + tar_name
        else:
            fname = self.dir + '/' + tar_name
                        
        if not tarfile.is_tarfile(fname):
            raise IOError(f'{fname} is not a valid tar.gz file. '
                          'Try again with `force_download=True`')
        with tarfile.open(fname, 'r:gz') as tar:
            print(f'Extracting {fname}...')
            tar.extractall(path = self.dir)

        if delete_tar and fname.exists():
            print("Deleting", fname)
            fname.unlink() 
            
            
    def data_from_model(self, data: tuple, 
                        time_now: Quantity["time"]) -> tuple:
        """
        Extract the stellar properties necessary for feedback
        from a BPASS stellar model file.
        
        Args:
            data                (tuple): Matrix loaded from the stellar model file
            time_now (Quantity["time"]): Current simulation time
            
        Returns:
            Feedback properties
        """
        # Replace NaNs by 0s
        data = np.nan_to_num(data)
        # We want to extract the following values
        time = data[:, 1] * u.yr
        ind_now = np.where((time_now - time) > 0 * u.yr)[-1][0]

        # Need to use column 36 for 50-0.9-1
        M1 = data[ind_now, 36] * u.Msun # 5 recommended in manual
        L1 = 10**data[ind_now, 4] * u.Lsun
        T1 = 10**data[ind_now, 3] * u.K
        R1 = 10**data[ind_now, 2] * u.Rsun

        M2 = data[ind_now, 37] * u.Msun
        L2 = 10**data[ind_now, 48] * u.Lsun
        T2 = 10**data[ind_now, 47] * u.K
        R2 = 10**data[ind_now, 46] * u.Rsun

        dM  = (data[ind_now, 39] + data[ind_now, 40] - data[ind_now, 41] - data[ind_now, 42] + \
               data[ind_now, 43] + data[ind_now, 44]) * u.Msun / (1.989 * u.s)
        
        return M1, M2, R1, R2, T1, T2, L1, L2, dM
    
    def read_bpass_data(self) -> tuple:
        """
        Read BPASS data from a stellar model file.
        """
        # General file name, to append
        gen_fname = self.dir + 'NEWBINMODS/NEWBINMODS/' + self.metal + \
                    '/sneplot-' + self.metal
        # Get frequency for each unique binary
        props = np.vstack((self.mass, self.mrat, self.logp))
        vals, _, counts = np.unique(props, return_index = True, 
                                    return_counts = True, axis=1)
        feedback = np.empty((10, len(counts)))
        feedback[0, :] = counts
        # Get the properties for each unique binary
        for i in range(len(counts)):
            _mass, _mrat, _logp = np.round(vals[:, i], decimals=1)
            # Matching file names from BPASS
            if str(_mass)[-2:] == '.0':
                _mass = int(_mass)
            if str(_logp)[-2:] == '.0':
                _logp = int(_logp)
            fname = gen_fname + '-' + str(_mass) + '-' + str(_mrat) + \
                    '-' + str(_logp)

            data = self.data_from_model(np.genfromtxt(fname), self.time)
            feedback[1, i] = data[0].to_value(u.Msun)
            feedback[2, i] = data[1].to_value(u.Msun)
            feedback[3, i] = data[2].to_value(u.Rsun)
            feedback[4, i] = data[3].to_value(u.Rsun)
            feedback[5, i] = data[4].to_value(u.K)
            feedback[6, i] = data[5].to_value(u.K)
            feedback[7, i] = data[6].to_value(u.Lsun)
            feedback[8, i] = data[7].to_value(u.Lsun)
            feedback[9, i] = data[8].to_value(u.Msun/u.yr)
             
        return feedback
            
            
    def get_stellar_models(self, tar_name = "bpass-v2.2-newmodels.tar.gz",
                           zip_name = "bpass_v2.2.zip", url = bpass_url) -> tuple :
        """
        Look for the stellar models for the stellar population
        and download/unzip/untar them if they are not available.
        """
        if os.path.isdir(self.dir + 'NEWBINMODS'):
            print("Reading BPASS data...")
            _tmp = self.read_bpass_data()
        else:
            print("BPASS data not available at", self.dir + 'NEWBINMODS')
            print("Looking for a tar file...")
            if os.path.isfile(self.dir + tar_name):
                self.untar(tar_name, delete_tar = False)
                _tmp = self.read_bpass_data()

            else:
                print("tar file not available at", self.dir)
                print("Looking for a zip file...")
                if os.path.isfile(self.dir + zip_name):
                    self.unzip(zip_name, target_file = tar_name, 
                               delete_zip = False, inspect = False)
                    self.untar(tar_name, delete_tar = False)
                    _tmp = self.read_bpass_data()

                else:
                    print("zip file not available at", self.dir)
                    if self.download:
                        self.downloader(url)
                        self.unzip(zip_name, target_file = tar_name, 
                               delete_zip = False, inspect = False)
                        self.untar(tar_name, delete_tar = False)
                        _tmp = self.read_bpass_data()

                    else:
                        print("Set force_download = True to download the files")
                        _tmp = np.empty((10, 1))


        return _tmp[0], _tmp[1]*u.Msun, _tmp[2]*u.Msun, _tmp[3]*u.Rsun, _tmp[4]*u.Rsun, \
               _tmp[5]*u.K, _tmp[6]*u.K, _tmp[7]*u.Lsun, _tmp[8]*u.Lsun, _tmp[9]*u.Msun / u.yr
