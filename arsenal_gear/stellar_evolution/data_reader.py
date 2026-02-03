"""
data_reader.py
==========

This file defines classes for reading stellar evolution data from
various sources into a uniform data structure format defined in
se_data_structures.py
"""

import time
from abc import ABC, abstractmethod
from pathlib import Path

import astropy.units as u
import numpy as np

from ..utils import downloader, extract_one, find_match, is_valid_txz
from .se_data_structures import Isochrone, IsochroneSet, StellarTrack, TrackSet


class IsochroneDataReader(ABC):
    """
    Abstract base class for reading isochrone data from various sources.
    """
    def __init__(self, **kwargs) -> None:
        # log10(Z/Zsun)
        self.met = kwargs.get('met', 0.0)
        # decides on verbose output
        self.verbose = kwargs.get('verbose', False)
        # decides whether or not to force a download of the isochrone data
        self.force_download = kwargs.get('force_download', False)

    @abstractmethod
    def read_iso_data(self) -> IsochroneSet:
        """
        Abstract method for reading isochrone data.
        """

    @abstractmethod
    def read_track_data(self) -> TrackSet:
        """
        Abstract method for reading stellar track data.
        """

class MISTReader(IsochroneDataReader):
    """
    Class for reading MIST isochrone and stellar track data.
    makes heavy use of code provided by the MIST team for reading their models.
    This is an instantiation of the IsochroneDataReader abstract base class.
    """
    # basic options for MIST isochrones
    mist_url = "https://waps.cfa.harvard.edu/MIST/data/tarballs_v1.2/{}"
    vcrits   = ["0.0", "0.4"]
    mets     = ["m4.00", "m3.50", "m3.00", "m2.50", "m2.00", "m1.75", "m1.50", "m1.25",\
                "m1.00", "m0.75", "m0.50", "m0.25", "p0.00", "p0.25", "p0.50"]
    def __init__(self, **kwargs) -> None:
        """
        Args:
            kwargs: Keyword arguments for the isochrone system.

        Usage:
            >> iso = MIST(met=0.0, vvcrit="0.0", interp_op="track")
            >> mini = np.array([1.0, 5.0, 10.0])*u.Msun
            >> lbol = iso.lbol(mini, 10*u.Myr)
            >> teff = iso.teff(mini, 10*u.Myr)

        Attributes:
            version     Dictionary containing the MIST and MESA version numbers.
            abun        Dictionary containing Yinit, Zinit, [Fe/H], and [a/Fe] values.
            rot         Rotation in units of surface v/v_crit.
            hdr_list    List of column headers.

        Methods:
            get_data        Retrieves the MIST data
            read_iso_data   Processes isochrone data into an IsochroneSet data structure
            read_track_data Processes stellar track data into a TrackSet data structure
        """
        # set input parameters
        super().__init__(**kwargs)
        self.vvcrit = kwargs.get("vvcrit", "0.0")
        self.rootdir = kwargs.get("rootdir", None)
        # Fill out dummy attributes
        self.hdr_list = None
        self.num_ages = None
        self.version = None
        self.rot = None
        self.abun = None

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

        # set directories, download data if necessary
        data_acq_start = time.time()
        self.get_data()
        data_acq_end = time.time()
        if self.verbose:
            print("Time to acquire data: ", data_acq_end - data_acq_start)

    def get_data(self):
        """
        Retrieves the MIST data from a local directory or downloads it if necessary.
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
        elif self.interp_op == "track":
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
            modeldir_path = self.rootdir / self.modeldir
            tarfile_path = self.rootdir / self.tarfile
            if not modeldir_path.is_dir():
                # model directory doesn't exist -> check for tarfile
                if not tarfile_path.is_file():
                    # if it doesn't exist, download it
                    force_download = True
                elif not is_valid_txz(tarfile_path):
                    # if it's not a valid tar file, download it again
                    force_download = True
                else:
                    # tarfile exists and is valid -> extract it
                    extract_one(tarfile_path, self.rootdir, delete_txz=True)
            elif ((self.interp_op == "iso") and
                  (not (modeldir_path / self.isofile).is_file())):
                # model directory exists but isochrone file doesn't -> force download
                force_download = True
        self.force_download = force_download
        if self.verbose:
            print("Force download: ", self.force_download)

        if self.force_download:
            url = self.mist_url.format(self.tarfile)
            message = "Downloading MIST data"
            # download the tarfile
            downloader(self.rootdir / self.tarfile, url, message)
            extract_one(self.rootdir / self.tarfile, self.rootdir, delete_txz=True)

    def read_iso_data(self) -> IsochroneSet:
        """
        Reads in isochrone data for MIST
        """
        lages, hdr_list, isos_mist = self.read_iso_file()
        isos = []
        for iso in isos_mist:
            age = np.power(10, iso["log10_isochrone_age_yr"][0])*u.yr
            qs = {}
            for label in hdr_list:
                if label != "log10_isochrone_age_yr":
                    qs[label] = iso[label]
            iso_data = Isochrone(age=age,
                                    eep_name='EEP',
                                    mini_name='initial_mass',
                                    lteff_name='log_Teff',
                                    llbol_name='log_L',
                                    qs=qs)
            isos.append(iso_data)
        self.hdr_list = hdr_list
        # get the maximum mass still alive for each isochrone
        max_mass = np.max(isos_mist[0]['initial_mass'])
        iset = IsochroneSet(lages=lages,
                            hdr_list=self.hdr_list,
                            isos=isos,
                            max_mass=max_mass*u.Msun)
        return iset

    def read_track_data(self) -> TrackSet:
        """
        Reads in stellar track data for MIST
        """
        eep_file_pattern = [("?????M.track.eep",),]
        mass_file_list = find_match(eep_file_pattern, self.rootdir / self.modeldir)

        masses = []
        tracks = []
        min_ages = []
        max_ages = []
        max_eep = 0
        for f in mass_file_list:
            # read in the EEP file
            minit, eeps, min_age, max_age, data = self.read_eep_file(f)
            # store the mass and EEP data
            max_eep = max(max_eep, np.max(eeps))
            masses.append(minit)
            min_ages.append(min_age)
            max_ages.append(max_age)
            tracks.append(StellarTrack(mass=minit*u.Msun,
                                        eeps=eeps,
                                        age_name='star_age',
                                        lteff_name='log_Teff',
                                        llbol_name='log_L',
                                        qs=data))
        masses = np.array(masses)*u.Msun
        min_ages = np.array(min_ages)*u.yr
        max_ages = np.array(max_ages)*u.yr
        self.hdr_list = data.keys()
        tset = TrackSet(masses=masses,
                                min_ages=min_ages,
                                max_ages=max_ages,
                                hdr_list=self.hdr_list,
                                tracks=tracks,
                                max_eep=max_eep)
        return tset

    def read_iso_file(self):
        """
        Reads in the isochrone file.

        This is essentially copied from read_mist_models.py in the MIST Github
        https://github.com/jieunchoi/MIST_codes/blob/master/scripts/read_mist_models.py

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

        This is essentially copied from read_mist_models.py in the MIST Github
        https://github.com/jieunchoi/MIST_codes/blob/master/scripts/read_mist_models.py

        """

        evol = np.loadtxt(fname, skiprows=11).T
        with open(fname, encoding='utf-8') as f:
            content = [line.split() for line in f]

        if not hasattr(self, 'version'):
            self.version = {'MIST': content[0][-1], 'MESA': content[1][-1]}
            self.abun = {content[3][i]:float(content[4][i]) for i in range(1,5)}
            self.rot = float(content[4][-1])
            self.hdr_list = content[11][1:]
        minit = float(content[7][1])
        hdr_list = content[11][1:]
        eeps = [int(j) for j in content[8][2:]]


        data = {hdr_list[i]:evol[i] for i in range(len(hdr_list))}
        min_age = min(data["star_age"])
        max_age = max(data["star_age"])

        return minit, eeps, min_age, max_age, data
