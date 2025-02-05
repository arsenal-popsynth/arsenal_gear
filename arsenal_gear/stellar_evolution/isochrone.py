"""
isochrone
==========

This submodule defines the interface to various stellar evolution codes
through interpreting and processing their isochrones
"""

import os

import astropy.units as u
import numpy as np
from astropy.units import Quantity


class Isochrone():
    """
    This class is used to load and interpret isochrones from various sources
    """
    def __init__(self) -> None:
        # laying out the basic attributed of the isochrone class
        self.filename = ""
        # the metallicity of the isochrone relative to solar
        self.metallicity = 1.0


    def get_Mmax(self, age:Quantity["time"]) -> Quantity["mass"]:
        """
        get the maximum mass of the stellar population that hasn't
        died yet (in e.g. a SN) as a function of age
        """

class MIST(Isochrone):
    """

    Reads in MIST isochrone files.


    """
    def __init__(self, filename : str, verbose:bool=False) -> None:
        """
        Args:
            filename: the name of .iso file.

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
        self.filename = filename
        if verbose:
            print('Reading in: ' + self.filename)

        self.ages, self.hdr_list, self.isos = self.read_iso_file()
        self.metallicity = self.abun['[Fe/H]']
        super().__init__()

    def read_iso_file(self):
        """

        Reads in the isochrone file.

        Args:
            filename: the name of .iso file.

        """
        print(os.getcwd())

        #open file and read it in
        with open(self.filename, encoding='utf-8') as f:
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
