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
                 bpass_dir: str, force_download: bool = False) -> None:
        """
        Args:
            binaries:  the BinaryPopulation instance for which we want information
            metal:     the metallicity of the stellar population
            TO DO -CCC, 26/06/2025: Read metallicity from binaries
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
            
        self.mass:    list = binaries.primary["mass"].to(u.Msun).value
        self.mrat:    list = binaries.primary["mass"] / binaries.secondary["mass"]
        self.logp:    list = np.log10(binaries["period"].to(u.d).value)
        super().__init__()
                

    def downloader(self, url = bpass_url, message = "Trying to download stellar models..."):
        """
        Method for downloading BPASS stellar models from the web.

        Args:
            url (str): The URL of the BPASS repository.
            directory (Path): The path of the directory to save the downloaded file.
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
                
    def unzip(self, zip_name, target_file=None, delete_zip=False, inspect=False):
        """
        Un-compress a zip file.
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
                
    
    def untar(self, tar_name, delete_tar=False):
        """
        Un-compress the stellar models tar.gz file.
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
    

                
class BPASS_population():
    """
    Reads in BPASS stellar model files for pre-processed BPASS populations.
    """
    # basic options for BPASS stellar models
    bpass_url = "https://www.dropbox.com/scl/fo/mpuas1xh5owmdadu0vpev/h?dl=0&e=1&rlkey=7vlk7ra6kvoztzmae8wr34kmz"
    population_to_download = ""
    
    
    def __init__(self, binaries: BinaryPopulation,
                 model: str, metal: str, bpass_dir: str) -> None:
        """
        Args:
            binaries:  the BinaryPopulation instance for which we want information
            model:     the BPASS IMF model we want to use
            metal:     the metallicity of the stellar population
            TO DO -CCC, 26/06/2025: Read metallicity from binaries
            bpass_dir: the directory for the BPASS models

        """
        self.models = ['100_100', '100_300', '135_100', '135_300', '135all_100',
                       '170_100', '170_300', 'chab100', 'chab300']
        if model not in self.models:
            raise ValueError(f"Model {model} does not exist.")
            
        self.metals = ["zem5", "zem4", "z001", "z002", "z003", "z004", "z006",
                       "z008", "z010", "z014",  "z020",  "z030",  "z040"]
        if metal not in self.metals:
            raise ValueError(f"Metallicity {metal} does not exist.")
            
    @staticmethod
    def downloader(model, url, directory, message):
        """
        Method for downloading a BPASS population dataset from the web.

        Args:
            model (str or Path): The IMF model name.
            url (str): The URL of the BPASS repository.
            directory (Path): The path of the directory to save the downloaded file.
            message (str): Optional message to display before downloading.

        Raises:
            Exception: If the download fails.

        """
        if message is not None:
            print(message)

        download_url = url + 'bpass_v2.2.1_imf' + model + '.tar.gz'
        if directory[-1] == '/':
            fname = directory + 'bpass_v2.2.1_imf' + model + '.tar.gz'
        else:
            fname = directory + '/bpass_v2.2.1_imf' + model + '.tar.gz'
            
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
