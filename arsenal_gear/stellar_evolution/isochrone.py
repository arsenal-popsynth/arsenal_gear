"""
isochrone
==========

This submodule defines the interface to various stellar evolution codes
through interpreting and processing their isochrones
"""

from functools import reduce
import glob
import os
import os.path as osp
from pathlib import Path
import tarfile
import time

import requests
from tqdm import tqdm

import astropy.units as u
import numpy as np
from scipy.interpolate import CubicSpline, pchip_interpolate
from astropy.units import Quantity

from . import utils as se_utils

class Isochrone():
    """
    This class is used to load and interpret isochrones from various sources
    """
    def __init__(self, **kwargs) -> None:
        # log10(Z/Zsun)
        self.met = kwargs.get('met', 0.0)
        # determines whether or not thie isochrone instance
        # is being used for testing or not. This changes the 
        # selection of the isochrone data to leave out the 
        # data being compared against.
        self.test = kwargs.get('test', False)
        # decides on verbose output
        self.verbose = kwargs.get('verbose', False)
        # whether or not to print output related to code profiling
        self.profile = kwargs.get('profile', False)
        # decides whether or not to force a download of the isochrone data
        self.force_download = kwargs.get('force_download', False)

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

    @staticmethod
    def _find_match(patterns, base_dir=None):
        if base_dir is None:
            glob_match = lambda p: sorted(glob.glob(*p))
        else:
            glob_match = lambda p: sorted(glob.glob(osp.join(base_dir,*p)))
        for p in patterns:
            f = glob_match(p)
            if f:
                break

        return f

    def mmax(self, t:Quantity["time"]) -> Quantity["mass"]:
        """
        get the maximum mass of the stellar population that hasn't
        died yet (in e.g. a SN) as a function of age, t

        default function in base class
        """
        return 0.0*u.Msun

    def mmaxdot(self, t: Quantity["time"]) -> Quantity["mass"]:
        """
        get the rate of change of the maximum mass of the stellar population
        with respect to time.

        Default function in base class
        Args:
            t: the age of the isochrone.
        """
        return 0.0*u.Msun/u.Myr

    def lbol(self, mini:Quantity["mass"], t: Quantity["time"]) -> Quantity["power"]:
        """
        get the bolometric luminosity of a star of initial mass mini at age age

        Default function in base class
        
        Args:
            mini: the initial mass of the star.
            t: the age of the isochrone.
        Returns:
            Quantity["power"]: the bolometric luminosity of the star.
        """
        return 0.0*u.erg/u.s

class MIST(Isochrone):
    """

    Reads in MIST isochrone files.


    """
    # basic options for MIST isochrones
    mist_url = "https://waps.cfa.harvard.edu/MIST/data/tarballs_v1.2/{}"
    vcrits   = ["0.0", "0.4"]
    mets     = ["m4.00", "m3.50", "m3.00", "m2.50", "m2.00", "m1.75", "m1.50", "m1.25",\
                "m1.00", "m0.75", "m0.50", "m0.25", "p0.00", "p0.25", "p0.50"]
    def __init__(self, **kwargs) -> None:
        """
        Args:
            met: the metallicity of the isochrone relative to solar

        Usage:
            >> iso = read_mist_models.ISO('MIST_v1.0_feh_p0.00_afe_p0.0_vvcrit0.4.iso')
            >> age_ind = iso._age_index(8.0)
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
        super().__init__(**kwargs)
        self.vvcrit = kwargs.get("vvcrit", "0.0")
        self.rootdir = kwargs.get("rootdir", None)

        if self.met<0:
            self.metstr = f"m{-1*self.met:.2f}"
        else:
            self.metstr = f"p{self.met:.2f}"
        if self.metstr not in self.mets:
            raise ValueError("Metallicity must be one of: " + str(self.mets))
        if self.vvcrit not in ["0.0", "0.4"]:
            raise ValueError("vvcrit must be either 0.0 or 0.4 for MIST isochrones")

        # decide whether to use EEPs or isochrones to interpolate
        # creating between for new isochrones
        self.interp_op = kwargs.get("interp_op", "iso")
        if self.interp_op not in ["iso", "eep"]:
            raise ValueError("interp_op must be either iso or eep")

        # set directories, download data if necessary
        data_acq_start = time.time()
        self.get_data()
        data_acq_end = time.time()
        if self.verbose:
            print("Time to acquire data: ", data_acq_end - data_acq_start)

        if self.verbose:
            print('Reading in data...')
        data_read_start = time.time()
        self.read_data()
        data_read_end = time.time()
        if self.verbose:
            print("Time to read data: ", data_read_end - data_read_start)

    def get_data(self):
        """
        Retrieves the MIST isochrone data from a local directory or downloads it if
        necessary.

        Returns:
            None
        """
        # set model directory, tarfile, and filename
        if self.interp_op == "iso":
            mod_temp = "MIST_v1.2_vvcrit{}_basic_isos"
            self.modeldir = mod_temp.format(self.vvcrit)
            tar_temp = "MIST_v1.2_vvcrit{}_basic_isos.txz"
            self.tarfile = tar_temp.format(self.vvcrit)
            iso_fname_temp = "MIST_v1.2_feh_{}_afe_p0.0_vvcrit{}_basic.iso"
            self.isofile = iso_fname_temp.format(self.metstr, self.vvcrit)
        elif (self.interp_op == "eep"):
            mod_temp = "MIST_v1.2_feh_{}_afe_p0.0_vvcrit{}_EEPS"
            self.modeldir = mod_temp.format(self.metstr, self.vvcrit)
            tar_temp = "MIST_v1.2_feh_{}_afe_p0.0_vvcrit{}_EEPS.txz"
            self.tarfile = tar_temp.format(self.metstr, self.vvcrit)
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
                else:
                    # tarfile exists and is valid -> extract it
                    self.extract_one(self.rootdir / self.tarfile, self.rootdir, delete_txz=True)
            elif (self.interp_op == "iso") and (not (self.rootdir / self.modeldir / self.isofile).is_file()):
                # model directory exists but isochrone file doesn't -> force download
                force_download = True
        self.force_download = force_download
        if self.verbose:
            print("Force download: ", self.force_download)

        if self.force_download:
            url = self.mist_url.format(self.tarfile)
            message = "Downloading MIST data"
            # download the tarfile
            self.downloader(self.rootdir / self.tarfile, url, message)
            self.extract_one(self.rootdir / self.tarfile, self.rootdir, delete_txz=True)

    def read_data(self):
        if self.interp_op == "iso":
            self.ages, self.hdr_list, self.isos = self.read_iso_file()
            # get the maximum mass still alive for each isochrone
            na = self.num_ages
            self.mmaxes = np.array([np.max(self.isos[i]['initial_mass']) for i in range(na)])
            self.MMAX = self.mmaxes[0]
            self.metallicity = self.abun['[Fe/H]']
        else:
            eep_file_pattern = [("?????M.track.eep",),]
            mass_file_list = self._find_match(eep_file_pattern, self.rootdir / self.modeldir)
            mass_nums = [int(f.split("/")[-1].split("M")[0]) for f in mass_file_list]
            self.num_masses = len(mass_nums)

            masses = []
            eeps_list = []
            mass_set = []
            min_ages = []
            max_ages = []
            for f in mass_file_list:
                # read in the EEP file
                minit, eeps, min_age, max_age, data = self.read_eep_file(f)
                # store the mass and EEP data
                masses.append(minit)
                eeps_list.append(eeps)
                min_ages.append(min_age)
                max_ages.append(max_age)
                mass_set.append(data)
            self.masses = np.array(masses)
            self.eeps_list = eeps_list
            self.min_ages = np.array(min_ages)
            self.max_ages = np.array(max_ages)
            self.mass_set = mass_set
            self.MMAX = np.max(masses)

    def read_iso_file(self):
        """
        Reads in the isochrone file.

        """

        # open file and read it in
        fname = self.rootdir / self.modeldir / self.isofile
        with open(fname, encoding='utf-8') as f:
            content = [line.split() for line in f]
        if not hasattr(self, 'version'):
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
            # grab info for each isochrone
            num_eeps = int(data[counter][-2])
            num_cols = int(data[counter][-1])
            hdr_list = data[counter+2][1:]
            formats = tuple([np.int32]+[np.float64 for i in range(num_cols-1)])
            iso = np.zeros((num_eeps),{'names':tuple(hdr_list),'formats':tuple(formats)})
            # read through EEPs for each isochrone
            for eep in range(num_eeps):
                iso_chunk = data[3+counter+eep]
                iso[eep] = tuple(iso_chunk)
            iso_set.append(iso)
            ages.append(iso[0][1])
            counter+= 3+num_eeps+2
        return ages, hdr_list, iso_set

    def read_eep_file(self, fname):

        """
        Reads in an EEP file.

        """

        evol = np.loadtxt(fname, skiprows=11).T
        with open(fname) as f:
            content = [line.split() for line in f]

        if not hasattr(self, 'version'):
            self.version = {'MIST': content[0][-1], 'MESA': content[1][-1]}
            self.abun = {content[3][i]:float(content[4][i]) for i in range(1,5)}
            self.rot = float(content[4][-1])
            self.hdr_list = content[11][1:]
        minit = float(content[7][1])
        hdr_list = content[11][1:]
        num_eeps = int(content[7][2])
        eeps = [int(j) for j in content[8][2:]]


        data = {hdr_list[i]:evol[i] for i in range(len(hdr_list))}
        min_age = min(data["star_age"])
        max_age = max(data["star_age"])

        return minit, eeps, min_age, max_age, data

    def _age_index(self, age : Quantity["time"]) -> int:
        """

        Returns the index of the isochrone closest to the requested age
        that is also younger than the requested age.

        Args:
            age: the age of the isochrone.

        """
        lage = np.log10(age.to(u.yr).value)
        if (max(lage) > max(self.ages)) or (min(lage) < min(self.ages)):
            print('The requested age is outside the range. Try between '
                  + str(min(self.ages)) + ' and ' + str(max(self.ages)))
            raise ValueError("Age is outside the range of the isochrones")

        ais = [np.where(np.array(self.ages) - la < 0)[0][-1] for la in lage]
        ais = np.array(ais, dtype=int)

        return ais

    def _get_ai_range(self, ai:int, n:int)-> int:
        """
        Gets the range of age indices to use for isochrone interpolation
        on n nearby points around age index ai. Leaves out the nearest
        isochrone if in testing mode.
        Args:
            ai: the index of the isochrone.
            n: the number of nearby points to use for interpolation.
        Returns:
            ais: the indices of the isochrones to use for interpolation.
        """
        nages = len(self.ages)

        if (self.test and ((ai==0) or (ai == nages-1))):
            raise ValueError("Test Isochrone must be in middle of isochrone range")
        if (n<=1):
            raise ValueError("Must Request more than one point for interpolation")
    
        # decide on right range of isochrone indices
        if (ai-n//2 <= 0):
            if self.test:
                ais = np.concatenate((np.arange(0, ai), np.arange(ai+1, n+1)))
            else:
                ais = np.arange(0, n)
        elif (ai+n//2+1 >= nages):
            if self.test:
                ais = np.concatenate((np.arange(nages-n-1, ai), np.arange(ai+1, nages)))
            else:
                ais = np.arange(nages-n, nages)
        else:
            if n%2 == 0:
                if self.test:
                    ais = np.concatenate((np.arange(ai-n//2, ai),np.arange(ai+1, ai+n//2+1)))
                else:
                    ais = np.arange(ai-n//2+1, ai+n//2+1)
            else:
                if self.test:
                    ais = np.concatenate((np.arange(ai-n//2, ai),np.arange(ai+1,ai+n//2+2)))
                else:
                    ais = np.arange(ai-n//2, ai+n//2+1)
        return np.array(ais,dtype=int)

    @staticmethod
    def _fixed_eep_q(j:int, eeps:list, qs:list):
        """
        Returns the value of a isochrone quantity at a fixed EEP
        across severl isochrones at different times.
        Args:
            j: the index of the EEP to get values for
            eeps: the list of EEPs for each isochrone.
            qs: the list of quantities for each isochrone.
        Returns:
            qj: the quantity at the fixed EEP.
        """
        return [q[np.where(eep == j)[0][0]] for (q,eep) in zip(qs,eeps)]

    def _interp_iso_quantity(self, t:Quantity["time"], label:str,
                             method:str="pchip", make_monotonic:bool=False) -> (np.float64, np.float64):
        """
        Interpolates between isochrones using the EEP values to create an
        intermediate age isochrone
        Args:
            t: The age of the isochrone to be generated, should a single value.
            label: The label of the quantity to be interpolated.
        Returns:
            qi: The specified quantity at the requested age. This is a function of
                initial mass which can also be generated in this way.
            eepi: The EEP values corresponding to the interpolated isochrone.
        """
        if not(np.isscalar(t.value)):
            # throw error
            raise ValueError("t must be a scalar")
        else:
            t = np.array([t.value])*t.unit
    
        start = time.time()
        ai = self._age_index(t)[0]
        ais = self._get_ai_range(ai, 4)
        lages = [self.ages[i] for i in ais]
        lt = np.log10(t.to(u.yr).value)
        qs = [self.isos[i][label] for i in ais]
        eeps = [self.isos[i]['EEP'] for i in ais]

        # eeps present in all isochrones
        eepi = reduce(np.intersect1d, tuple(eeps))
        if (method == "pchip"):
            f = lambda x,x0,y0: pchip_interpolate(x0,y0,x)[0]
        elif (method == "linear"):
            f = lambda x,x0,y0: np.interp(x,x0,y0)[0]
        else:
            raise ValueError("method must be either pchip or linear")
        end = time.time()
        if self.profile:
            print("\t\tSet up of inerpolation took: ", end-start)

        # interpolate in log(age) at each eep
        start = time.time()
        qi = np.array([f(lt, lages, self._fixed_eep_q(j,eeps,qs)) for j in eepi])
        end = time.time()
        if self.profile:
            print("\t\tInterpolation took: ", end-start)


        if make_monotonic:
            start = time.time()
            if np.any(np.diff(qi) < 0):
                se_utils._make_monotonic(eepi,qi)
            end = time.time()
            if self.profile:
                print("\t\tMonotonic interpolation took: ", end-start)

        return (eepi,qi)

    def _get_mmax_age_interp(self):
        """
        Returns the maximum mass as a function of age,
        properly interpreted from either the isochrone or EEP data.
        """
        if self.interp_op == "iso":
            lages = self.ages
            mmax = self.mmaxes
        else:
            # interpreting from EEPs
            lages = np.log10(self.max_ages)[::-1]
            mmax = self.masses[::-1]
            sel = se_utils._index_monotonic(lages)
            lages = lages[sel]
            mmax = mmax[sel]
        return (lages, mmax)

    def mmax(self, t: Quantity["time"]) -> Quantity["mass"]:
        """
        get the maximum mass of the stellar population that hasn't
        died yet (in e.g. a SN) as a funciton of age, using a cubic spline
        based on maximum mass reported in the isochrone data
        """
        if np.isscalar(t.value):
            t = np.array([t.value])*t.unit

        lt = np.log10(t.to(u.yr).value)

        (lages, mmax) = self._get_mmax_age_interp()
        # cubic spline interpolation in log(age)
        interp = CubicSpline(lages,mmax,bc_type="clamped")(lt)
        interp[np.where(lt < min(lages))] = self.MMAX
        interp[np.where(lt > max(lages))] = 0.0
        return interp*u.Msun

    
    def mmaxdot(self, t: Quantity["time"]) -> Quantity["mass"]:
        """
        get the rate of change of the maximum mass of the stellar population
        with respect to time. Uses a cubic spline and takes the derivative
        """
        if np.isscalar(t.value):
            t = np.array([t.value])*t.unit

        lt = np.log10(t.to(u.yr).value)

        (lages, mmax) = self._get_mmax_age_interp()
        # return the first derivative of the cubic spline
        unitfac = u.Msun/t.to(u.Myr)/np.log(10)
        interp = CubicSpline(lages,mmax,bc_type="clamped")(lt, 1)*unitfac
        interp *= np.logical_and(lt > min(lages), lt < max(lages)) 
        return interp

    def lbol(self, mini:Quantity["mass"], t: Quantity["time"]) -> Quantity["power"]:
        """
        get the bolometric luminosity of a star of initial mass mini at age t
        Args:
            mini: the initial mass of the star. Can be an array
            t: the age of the isochrone. Should be a single time (for now...)
        Returns:
            Quantity["power"]: the bolometric luminosity of the star.
        """

        if np.isscalar(mini.value):
            mini = np.array([mini.value])*mini.unit
        
        # construct isochrone for mass/luminosity relationship
        start = time.time()
        (eepi,logLi) = self._interp_iso_quantity(t, 'log_L')
        end = time.time()
        if self.profile:
            print("\tTime to interpolate Lbol: ", end-start)
        start = time.time()
        (eepi,massi) = self._interp_iso_quantity(t, 'initial_mass', make_monotonic=True)
        end = time.time()
        if self.profile:
            print("\tTime to interpolate mass: ", end-start)

        logLbol_res = pchip_interpolate(massi, logLi, mini.value)
        Lbol_res = np.power(10, logLbol_res)*u.Lsun
        Lbol_res *= np.logical_and(mini.value > min(massi), mini.value < max(massi))
        return Lbol_res

    def teff(self, mini:Quantity["mass"], t: Quantity["time"],method:str="pchip") -> Quantity["temperature"]:
        """
        get the atmospheric effective temperature of a star of initial mass mini
        at age t
        Args:
            mini: the initial mass of the star. Can be an array
            t: the age of the isochrone. Should be a single time (for now...)
        Returns:
            Quantity["power"]: the bolometric luminosity of the star.
        """

        if np.isscalar(mini.value):
            mini = np.array([mini.value])*mini.unit

        # construct isochrone for mass/temperature relationship
        start = time.time()
        (eepi,logTi) = self._interp_iso_quantity(t, 'log_Teff',method=method)
        end = time.time()
        if self.profile:
            print("\tTime to interpolate Teff: ", end-start)
        start = time.time()
        (eepi,massi) = self._interp_iso_quantity(t, 'initial_mass',make_monotonic=True)
        end = time.time()
        if self.profile:
            print("\tTime to interpolate mass: ", end-start)

        logTeff_res = pchip_interpolate(massi, logTi, mini.value)
        Teff_res = np.power(10, logTeff_res)*u.K
        Teff_res *= np.logical_and(mini.value > min(massi), mini.value < max(massi))
        return Teff_res

    def mini(self, t: Quantity["time"],method:str="pchip") -> Quantity["temperature"]:
        """
        get the atmospheric effective temperature of a star of initial mass mini
        at age t
        Args:
            mini: the initial mass of the star. Can be an array
            t: the age of the isochrone. Should be a single time (for now...)
        Returns:
            Quantity["power"]: the bolometric luminosity of the star.
        """

        (eepi,massi) = self._interp_iso_quantity(t, 'initial_mass')
        return (eepi,massi)