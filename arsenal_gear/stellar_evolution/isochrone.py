"""
isochrone
==========

This submodule defines the interface to various stellar evolution codes
through interpreting and processing their isochrones
"""

import os
from pathlib import Path
import tarfile
import time

import requests
from tqdm import tqdm

import astropy.units as u
import numpy as np
from numpy.polynomial.polynomial import Polynomial
from scipy.interpolate import lagrange
from astropy.units import Quantity


class Isochrone():
    """
    This class is used to load and interpret isochrones from various sources
    """
    def __init__(self, met:np.float32) -> None:
        # laying out the basic attributed of the isochrone class
        self.filename = ""
        # the metallicity of the isochrone relative to solar
        self.met = met

    @staticmethod
    def downloader(fname, url, message):
        """
        Method for downloading isochrone data from the web where available.

        Args:
            fname (str or Path): The file name or path to save the downloaded file.
            url (str): The URL of the file to download.
            message (str): Optional message to display before downloading.

        Raises:
            Exception: If the download fails.

        """
        if message is not None:
            print(message)

        try:
            response = requests.get(url, stream=True, timeout=10)
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

    @staticmethod
    def is_valid_txz(fname):
        """
        Check if a .txz file is valid.

        Parameters:
        - fname (str): The path to the .txz file.

        Returns:
        - bool: True if the .txz file is valid, False otherwise.
        """
        try:
            with tarfile.open(fname, "r:xz") as tar:
                tar.getmembers()
            return True
        except (tarfile.TarError, ValueError, OSError, EOFError) as e:
            print(f"Invalid .txz file: {e}")
            return False

    def get_Mmax(self, age:Quantity["time"]) -> Quantity["mass"]:
        """
        get the maximum mass of the stellar population that hasn't
        died yet (in e.g. a SN) as a function of age
        """
        return 0.0*u.Msun

class MIST(Isochrone):
    """

    Reads in MIST isochrone files.


    """
    # basic options for MIST isochrones
    url_temp = "https://waps.cfa.harvard.edu/MIST/data/tarballs_v1.2/MIST_v1.2_vvcrit{}_basic_isos.txz"
    vcrits   = ["0.0", "0.4"]
    mets     = ["m4.00", "m3.50", "m3.00", "m2.50", "m2.00", "m1.75", "m1.50", "m1.25",\
                "m1.00", "m0.75", "m0.50", "m0.25", "p0.00", "p0.25", "p0.50"]
    def __init__(self, met:np.float32, vvcrit:str="0.0", rootdir:str=None,
                 force_download:bool=False, verbose:bool=False) -> None:
        """
        Args:
            met: the metallicity of the isochrone relative to solar

        Usage:
            >> iso = read_mist_models.ISO('MIST_v1.0_feh_p0.00_afe_p0.0_vvcrit0.4.iso')
            >> age_ind = iso.age_index(8.0)
            >> logTeff = iso.isos[age_ind]['log_Teff']
            >> logL = iso.isos[age_ind]['log_L']
            >> plt.plot(logTeff, logL) #plot the HR diagram for logage = 8.0

        Attributes:
            version     Dictionary containing the MIST and MESA version numbers.
            abun        Dictionary containing Yinit, Zinit, [Fe/H], and [a/Fe] values.
            rot         Rotation in units of surface v/v_crit.
            ages        List of ages.
            num_ages    Number of isochrones.
            hdr_list    List of column headers.
            isos        Data.
        """
        # set input parameters
        super().__init__(met)
        if met<0:
            self.metstr = f"m{-1*met:.2f}"
        else:
            self.metstr = f"p{met:.2f}"
        if self.metstr not in self.mets:
            raise ValueError("Metallicity must be one of: " + str(self.mets))
        if vvcrit not in ["0.0", "0.4"]:
            raise ValueError("vvcrit must be either 0.0 or 0.4 for MIST isochrones")
        self.vvcrit = vvcrit
        self.rootdir = rootdir
        self.force_download = force_download
        self.verbose = verbose

        # set directories, download data if necessary
        data_acq_start = time.time()
        self.get_data()
        data_acq_end = time.time()
        if verbose:
            print("Time to acquire data: ", data_acq_end - data_acq_start)

        if verbose:
            print('Reading in: ' + self.isofile)

        self.ages, self.hdr_list, self.isos = self.read_iso_file()
        # get the maximum mass still alive for each isochrone
        ai = np.arange(len(self.ages))
        self.mmaxes = np.array([np.max(self.isos[i]['initial_mass']) for i in ai])
        self.metallicity = self.abun['[Fe/H]']

    def get_data(self):
        """
        Retrieves the MIST isochrone data from a local directory or downloads it if
        necessary.

        Returns:
            None
        """
        # set model directory, tarfile, and filename
        mod_temp = "MIST_v1.2_vvcrit{}_basic_isos"
        self.modeldir = mod_temp.format(self.vvcrit)
        tar_temp = "MIST_v1.2_vvcrit{}_basic_isos.txz"
        self.tarfile = tar_temp.format(self.vvcrit)
        iso_fname_temp = "MIST_v1.2_feh_{}_afe_p0.0_vvcrit{}_basic.iso"
        self.isofile = iso_fname_temp.format(self.metstr, self.vvcrit)

        # make data storage directory if it doesn't exist
        rootdir = self.rootdir
        if rootdir is None:
            rootdir = Path(__file__).parent.absolute() / 'data/mist'
            if not rootdir.is_dir():
                rootdir.mkdir(parents=True, exist_ok=True)
        self.rootdir = Path(rootdir)

        # check if the file exists and download otherwise
        force_download = self.force_download
        if not force_download:
            if not (self.rootdir / self.modeldir).is_dir():
                # model directory doesn't exist -> check for tarfile
                if not (self.rootdir / self.tarfile).is_file():
                    # if it doesn't exist, download it
                    force_download = True                
                elif not self.is_valid_txz(self.rootdir / self.tarfile):
                    # if it's not a valid tar file, download it again
                    force_download = True
            elif not (self.rootdir / self.modeldir / self.isofile).is_file():
                # model directory exists but isochrone file doesn't -> force download
                force_download = True
        self.force_download = force_download
        if self.verbose:
            print("Force download: ", self.force_download)

        if self.force_download:
            url = self.url_temp.format(self.vvcrit)
            message = "Downloading MIST isochrone data"
            # download the tarfile
            self.downloader(self.rootdir / self.tarfile, url, message)
            self.extract_one(self.rootdir / self.tarfile, self.rootdir, delete_txz=False)
    
    @staticmethod
    def extract_one(fname, extractdir, delete_txz=False):
        """
        Unzips a single ZIP file.
        """
        # Ensure output directory exists
        os.makedirs(extractdir, exist_ok=True)
        if not tarfile.is_tarfile(fname):
            raise IOError(f'{fname} is not a valid txz file. '
                          'Try again with `force_download=True`')
        with tarfile.open(fname, 'r:xz') as tar:
            tar.extractall(path=extractdir)

        if delete_txz and fname.exists():
            print(fname)
            fname.unlink()

        return extractdir

    def read_iso_file(self):
        """

        Reads in the isochrone file.

        Args:
            filename: the name of .iso file.

        """

        # open file and read it in
        fname = self.rootdir / self.modeldir / self.isofile
        with open(fname, encoding='utf-8') as f:
            content = [line.split() for line in f]
        self.version = {'MIST': content[0][-1], 'MESA': content[1][-1]}
        self.abun = {content[3][i]:float(content[4][i]) for i in range(1,5)}
        self.rot = float(content[4][-1])
        self.num_ages = int(content[6][-1])

        #read one block for each isochrone
        iso_set = []
        ages = []
        counter = 0
        data = content[8:]
        for _ in range(self.num_ages):
            #grab info for each isochrone
            num_eeps = int(data[counter][-2])
            num_cols = int(data[counter][-1])
            hdr_list = data[counter+2][1:]
            formats = tuple([np.int32]+[np.float64 for i in range(num_cols-1)])
            iso = np.zeros((num_eeps),{'names':tuple(hdr_list),'formats':tuple(formats)})
            #read through EEPs for each isochrone
            for eep in range(num_eeps):
                iso_chunk = data[3+counter+eep]
                iso[eep]=tuple(iso_chunk)
            iso_set.append(iso)
            ages.append(iso[0][1])
            counter+= 3+num_eeps+2
        return ages, hdr_list, iso_set

    def age_index(self, age : Quantity["time"]) -> int:
        """

        Returns the index of the isochrone closest to the requested age
        that is also younger than the requested age.

        Args:
            age: the age of the isochrone.

        """
        lage = np.log10(age.to(u.yr).value)
        ais = [np.where(np.array(self.ages) - la < 0)[0][-1] for la in lage]
        ais = np.array(ais, dtype=int)

        if (max(lage) > max(self.ages)) or (min(lage) < min(self.ages)):
            print('The requested age is outside the range. Try between '
                  + str(min(self.ages)) + ' and ' + str(max(self.ages)))

        return ais

    def get_Mmax(self, age: Quantity["time"]) -> Quantity["mass"]:
        """
        get the maximum mass of the stellar population that hasn't
        died yet (in e.g. a SN) as a funciton of age
        """
        if np.isscalar(age.value):
            age = np.array([age.value])*age.unit

        # linear interpolation in log(age) to get the maximum mass
        ai = self.age_index(age)
        lage = np.log10(age.to(u.yr).value)
        a_lo = np.array([self.ages[i] for i in ai])
        a_hi = np.array([self.ages[i+1] for i in ai])
        mmax_lo = np.array([np.max(self.isos[i]['initial_mass']) for i in ai])
        mmax_hi = np.array([np.max(self.isos[i+1]['initial_mass']) for i in ai])
        s = (mmax_hi - mmax_lo)/(a_hi - a_lo)
        o_s = mmax_lo - s*a_lo
        if len(age) == 1:
            return (s*lage + o_s)[0]*u.Msun
        return (s*lage + o_s)*u.Msun

    def get_Mmax2(self, age: Quantity["time"]) -> Quantity["mass"]:
        """
        get the maximum mass of the stellar population that hasn't
        died yet (in e.g. a SN) as a funciton of age
        """
        if np.isscalar(age.value):
            age = np.array([age.value])*age.unit

        # linear interpolation in log(age) to get the maximum mass
        ai = self.age_index(age)
        lage = np.log10(age.to(u.yr).value)
        a_lo = np.array([self.ages[i] for i in ai])
        a_hi = np.array([self.ages[i+1] for i in ai])
        mmax_lo = np.array([np.max(self.isos[i]['initial_mass']) for i in ai])
        mmax_hi = np.array([np.max(self.isos[i+1]['initial_mass']) for i in ai])
        s = (np.log10(mmax_hi) - np.log10(mmax_lo))/(a_hi - a_lo)
        o_s = np.log10(mmax_lo) - s*a_lo
        if len(age) == 1:
            return (10**(s*lage + o_s))[0]*u.Msun
        return (10**(s*lage + o_s))*u.Msun

    def get_Mmax3(self, age: Quantity["time"]) -> Quantity["mass"]:
        """
        get the maximum mass of the stellar population that hasn't
        died yet (in e.g. a SN) as a funciton of age
        """
        if np.isscalar(age.value):
            age = np.array([age.value])*age.unit

        # quartic interpolation in log(age) to get the maximum mass
        lages = self.ages
        ai = self.age_index(age)
        lage = np.log10(age.to(u.yr).value)
        mmax = self.mmaxes
        aranges = [lages[max(i-1,0):i+3] for i in ai]
        mranges = [mmax[max(i-1,0):i+3] for i in ai]
        fs = [Polynomial(lagrange(aranges[j],mranges[j]).coef[::-1]) for j in range(len(ai))]
        return np.array([fs[j](lage[j]) for j in range(len(ai))])*u.Msun