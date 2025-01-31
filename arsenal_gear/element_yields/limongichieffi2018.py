import astropy.units as u
from astropy.units import Quantity
import numpy as np

import os
import re
from typing import Type

from .yields import Yields, Source

class LimongiChieffi2018(Yields):
    """
    Include yields majority of elements (and isotopes) for massive stars 
    [13, 15, 20, 25, 30, 40, 60, 80, 120] Msun, with initial iron abundance 
    ([Fe/H]) [0, -1, -2, -3], and three different models for stellar 
    rotation [0, 150, 300] km/s.

    Download yields and place/link to directory LimongiChieffi2018 placed
    in the path of this file, e.g.,
    ln -s /path/to/downloaded/yields /path/to/arsenal_gear/element_yields/LimongiChieffi2018

    Yields available at http://orfeo.iaps.inaf.it/

    Reference: Limongi M & Chieffi A, 2018, ApJS, 237, 13L
    """

    def __init__(self, model: str ='R'):

        self.models = ['F', 'I', 'M', 'R']
        self.mass = np.array(
            [13.0, 15.0, 20.0, 25.0, 30.0, 40.0, 60.0, 80.0, 120.0])
        self.metal = np.array([3.236e-5, 3.2363e-4, 3.236e-3, 1.345e-2])
        self.feh = np.array([-3, -2, -1, 0])
        self.rot = np.array([0, 150, 300])
        self.ccsn_mmax = 25

        if model not in self.models:
            raise ValueError("Model does not exist.")

        self.filedir = os.path.dirname(os.path.realpath(__file__))
        
        self.yield_tablefile = self.filedir + \
            f'/LimongiChieffi2018/tab_{model}/tab_yieldstot_ele_exp.dec'
        self.wind_tablefile = self.filedir+f'/LimongiChieffi2018/tab_{model}/tabwind.dec'

        self.elements, self.atomic_num = self.get_element_list()

        # Stellar wind yields
        wind_yld = self.load_wind_yields()
        self.wind = Source(
            self.elements, [self.rot, self.metal, self.mass], wind_yld)

        # Core-collapse SNe yields
        ccsn_yld = self.load_ccsn_yields()
        self.ccsn = Source(
            self.elements, [self.rot, self.metal, self.mass], ccsn_yld)

    def ccsn_yields(self, elements, mass, metal, rot, interpolate='nearest', extrapolate=False):
        """ Returns yields [Msol] for elements from core-collapse supernovae given mass [Msol], 
            metallicity [mass fraction], and rotation [km/s].
        """
        return self.ccsn.get_yld(elements, [rot, metal, mass], interpolate=interpolate, extrapolate=extrapolate)

    def wind_yields(self, elements, mass, metal, rot, interpolate='nearest', extrapolate=False):
        """ Returns yields [Msol] for elements from pre-supernovae wind given mass [Msol], 
            metallicity [mass fraction], and rotation [km/s].
        """
        return self.wind.get_yld(elements, [rot, metal, mass], interpolate=interpolate, extrapolate=extrapolate)

    def yields(self, elements, mass, metal, rot, interpolate='nearest', extrapolate=False):
        """ Returning yields [Msol] for each element in elements, given
            stellar parameters mass [Msol], metal [mass fraction], and rot [km/s]."""
        args = (elements, mass, metal, rot, interpolate, extrapolate)
        return self.wind_yields(*args) + self.ccsn_yields(*args)

    def ccsn_mloss(self, mass, metal, rot, interpolate='nearest', extrapolate=False):
        """ Returning mass loss [Msol] as sum of all elements from core-collapse supernovae given mass [Msol], 
            metallicity [mass fraction], and rotation [km/s].
        """
        return self.ccsn.get_mloss([rot, metal, mass], interpolate=interpolate, extrapolate=extrapolate)

    def wind_mloss(self, mass, metal, rot, interpolate='nearest', extrapolate=False):
        """ Returning mass loss [Msol] as sum of all elements from pre-supernovae wind given mass [Msol], 
            metallicity [mass fraction], and rotation [km/s].
        """
        return self.wind.get_mloss([rot, metal, mass], interpolate=interpolate, extrapolate=extrapolate)

    def mloss(self, mass, metal, rot, source='all', interpolate='nearest', extrapolate=False):
        """ Returning mass loss [Msol] as sum of all elements, given
            stellar parameters mass [Msol], metal [mass fraction], and rot [km/s]."""

        return self.wind_mloss(mass, metal, rot,
                               source=source,
                               interpolate=interpolate,
                               extrapolate=extrapolate) \
            + self.ccsn_mloss(mass, metal, rot,
                              source=source,
                              interpolate=interpolate,
                              extrapolate=extrapolate)

    def get_element_list(self):

        with open(self.wind_tablefile, 'r') as file:
            lines = file.readlines()

        elements = []
        atomic_num = []
        for line in lines[1:]:
            if line.split()[0] == 'ele':
                return elements, np.array(atomic_num, dtype=float)
            element = ''.join(re.findall(r'[a-zA-Z]', line.split()[0]))
            if element not in elements:
                elements.append(element)
                atomic_num.append(int(line.split()[1]))

    def load_wind_yields(self):

        wind_yld = np.zeros([len(self.elements), self.rot.size,
                             self.metal.size, self.mass.size])

        with open(self.wind_tablefile, 'r') as file:
            lines = file.readlines()

        for index, line in enumerate(lines):
            if line.split()[0] == 'ele':
                model = line.split()[4]
                ind_metal = self.get_metal_index_from_model(model[3])
                ind_rot = self.get_rot_index_from_model(model[4:])

                data = np.genfromtxt(self.wind_tablefile,
                                     usecols=[1, 4, 5, 6, 7, 8, 9, 10, 11, 12],
                                     skip_header=index+1,
                                     max_rows=142).T
                for i, (atom_nr, element) in enumerate(zip(self.atomic_num, self.elements)):
                    mask = (data[0] == atom_nr)
                    wind_yld[i, ind_rot, ind_metal] = np.sum(
                        data[1:, mask], axis=1)

        return wind_yld

    def load_ccsn_yields(self):

        wind_yld = self.load_wind_yields()
        total_yld = np.zeros([len(self.elements), self.rot.size,
                              self.metal.size, self.mass.size])

        with open(self.yield_tablefile, 'r') as file:
            lines = file.readlines()

        for index, line in enumerate(lines):
            if line.split()[0] == 'ele':
                model = line.split()[4]
                ind_metal = self.get_metal_index_from_model(model[3])
                ind_rot = self.get_rot_index_from_model(model[4:])

                total_yld[:, ind_rot, ind_metal, :] = np.genfromtxt(self.yield_tablefile,
                                                                    usecols=[
                                                                        4, 5, 6, 7, 8, 9, 10, 11, 12],
                                                                    skip_header=index+1,
                                                                    max_rows=53)

        ccsn_yld = total_yld - wind_yld
        ccsn_yld[ccsn_yld < 0.0] = 0.0
        ccsn_yld[:, :, :, (self.mass > self.ccsn_mmax)] = 0.0

        return ccsn_yld

    @staticmethod
    def get_metal_index_from_model(model):
        if model == 'a':
            return 3
        elif model == 'b':
            return 2
        elif model == 'c':
            return 1
        elif model == 'd':
            return 0
        else:
            raise ValueError("Model does not exist.")

    @staticmethod
    def get_rot_index_from_model(model):
        if model == '000':
            return 0
        elif model == '150':
            return 1
        elif model == '300':
            return 2
        else:
            raise ValueError("Model does not exist.")
