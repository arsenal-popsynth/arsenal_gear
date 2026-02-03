"""
isochrone
==========

This submodule defines the interface to various stellar evolution codes
through interpreting and processing their isochrones
"""

from functools import reduce

import astropy.units as u
import numpy as np
from astropy.units import Quantity
from scipy.interpolate import pchip_interpolate

# local imports
from ..utils import array_utils, masked_power
from .se_data_structures import Isochrone
from .data_reader import MISTReader

class IsochroneInterpolator():
    """
    This class is used to interpolate isochrones from various sources
    It can be used to interpolate stellar tracks or isochrones from any 
    stellar evolution code, provided there is a IsochroneDataReader class
    implemented in order to read the data.
    """

    def __init__(self, **kwargs) -> None:
        # log10(Z/Zsun)
        self.met = kwargs.get('met', 0.0)
        # determines whether or not the isochrone instance
        # is being used for testing or not. This changes the 
        # selection of the isochrone data to leave out the 
        # data being compared against.
        self.test = kwargs.get('test', False)
        # decides whether or not to force a download of the isochrone data
        self.interp_op = kwargs.get("interp_op", "iso")
        if self.interp_op not in ["iso", "eep"]:
            raise ValueError("interp_op must be either iso or eep")
    
        self.isochrone_opt = kwargs.get("isochrone_opt", "mist")
        if self.isochrone_opt.lower() == "mist":
            self.reader = MISTReader(**kwargs)
        else:
            raise ValueError("isochrone_opt must be mist for now")

        if self.interp_op == "iso":
            self.iset = self.reader.read_iso_data()
            self.llbol_label = self.iset.isos[0].llbol_name
            self.lteff_label = self.iset.isos[0].lteff_name
        else:
            self.tset = self.reader.read_track_data()
            self.llbol_label = self.tset.tracks[0].llbol_name
            self.lteff_label = self.tset.tracks[0].lteff_name

    @staticmethod
    def _get_interpolator(method:str):
        """
        Returns the appropriate interpolation function based on the method string
        Args:
            method: the interpolation method to use, either pchip or linear
        Returns:
            interp: the interpolation function
        """
        if method == "pchip":
            def interp(x,x0,y0):
                return pchip_interpolate(x0,y0,x)
        elif method == "linear":
            def interp(x,x0,y0):
                return np.interp(x,x0,y0)
        else:
            raise ValueError("method must be either pchip or linear")
        return interp

    def supplement_labels(self, labels:list[str]) -> list[str]:
        """
        Adds basic labels to the list of isochrone quantities to be interpolated.
        Args:
            labels: The list of labels to be interpolated.
        Returns:
            labels: The supplemented list of labels including basic quantities.
        """
        if self.interp_op == "iso":
            necessary_labels = [self.iset.isos[0].mini_name,
                                self.iset.isos[0].lteff_name,
                                self.iset.isos[0].llbol_name]
        else:
            necessary_labels = [self.tset.tracks[0].lteff_name,
                                self.tset.tracks[0].llbol_name]
        for nl in necessary_labels:
            if nl not in labels:
                labels.append(nl)
        return labels

    ######## ISOCHRONE INTERPOLATION FUNCTIONS ########
    def _age_index(self, age: Quantity["time"]) -> int:
        """

        Returns the index of the isochrone closest to the requested age
        that is also younger than the requested age.

        Args:
            age: the age of the isochrone.
        Returns:
            ai: the index of the isochrone in the isochrone set
        """
        lage = np.log10(age.to(u.yr).value)
        ages = self.iset.lages
        if (max(lage) >= max(ages)) or (min(lage) <= min(ages)):
            print('The requested age is outside the range. Try between '
                  + str(min(ages)) + ' and ' + str(max(ages)))
            raise ValueError("Age is outside the range of the isochrones")

        ais = [np.where(np.array(ages) - la < 0)[0][-1] for la in lage]
        ais = np.array(ais, dtype=int)

        return ais

    def _get_ai_range(self, ai: int, n: int) -> int:
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
        nages = len(self.iset.lages)

        if self.test and (ai in (0, nages - 1)):
            raise ValueError("Test Isochrone must be in middle of isochrone range")
        if n <= 1:
            raise ValueError("Must Request more than one point for interpolation")

        # decide on right range of isochrone indices
        if ai - n // 2 <= 0:
            if self.test:
                ais = np.concatenate((np.arange(0, ai), np.arange(ai + 1, n + 1)))
            else:
                ais = np.arange(0, n)
        elif ai + n // 2 + 1 >= nages:
            if self.test:
                ais = np.concatenate(
                    (np.arange(nages - n - 1, ai), np.arange(ai + 1, nages))
                )
            else:
                ais = np.arange(nages - n, nages)
        else:
            if n % 2 == 0:
                if self.test:
                    ais = np.concatenate(
                        (np.arange(ai - n // 2, ai), np.arange(ai + 1, ai + n // 2 + 1))
                    )
                else:
                    ais = np.arange(ai - n // 2 + 1, ai + n // 2 + 1)
            else:
                if self.test:
                    ais = np.concatenate(
                        (np.arange(ai - n // 2, ai), np.arange(ai + 1, ai + n // 2 + 2))
                    )
                else:
                    ais = np.arange(ai - n // 2, ai + n // 2 + 1)
        return np.array(ais, dtype=int)

    @staticmethod
    def _fixed_eep_q(j: int, eeps: list, qs: list):
        """
        Returns the value of a isochrone quantity at a fixed EEP
        across several isochrones at different times.
        Args:
            j: the index of the EEP to get values for
            eeps: the list of EEPs for each isochrone.
            qs: the list of quantities for each isochrone.
        Returns:
            qj: the quantity at the fixed EEP.
        """
        return [q[np.where(eep == j)[0]][0] for (q, eep) in zip(qs, eeps)]

    def _construct_iso_isochrone(self, t:Quantity["time"], labels:list[str],
                                 method:str="pchip",
                                 supplement=True) -> Isochrone:
        """
        Interpolates between isochrones using the EEP values to create an
        intermediate age isochrone
        Args:
            t: The age of the isochrone to be generated, should be a single value.
            labels: The labels of the quantities to be interpolated.
            method: the interpolation method to use, either pchip or linear
            supplement: Whether to supplement the labels with basic quantities
        Returns:
            iso: the isochrone at the requested age t with the requested quantities
                 given by the labels list. 
        """
        if supplement:
            labels = self.supplement_labels(labels)

        if np.isscalar(t.value):
            ai = self._age_index(Quantity([t.value], t.unit))[0]
        else:
            if len(t) != 1:
                raise ValueError("t must be a single value for now")
            ai = self._age_index(Quantity([t.value[0]], t.unit))[0]
        # get nearby isochrones and interpolate between them
        ais = self._get_ai_range(ai, 4)
        lages = np.array([self.iset.lages[i] for i in ais])
        lt = np.log10(t.to(u.yr).value)
        qs = [[self.iset.isos[i].qs[label] for i in ais] for label in labels]
        eep_str = self.iset.isos[0].eep_name
        eeps = [self.iset.isos[i].qs[eep_str] for i in ais]

        # eeps present in all isochrones
        eepi = reduce(np.intersect1d, tuple(eeps))
        f = self._get_interpolator(method)

        # interpolate in log(age) at each eep
        qis = [np.array([f(lt, lages, self._fixed_eep_q(j,eeps,q)) for j in eepi]) for q in qs]

        # make sure initial mass is monotonic in interpolation
        for (i,label) in enumerate(labels):
            if label == self.iset.isos[0].mini_name:
                if np.any(np.diff(qis[i]) <= 0):
                    qis[i] = array_utils.make_monotonic_increasing(eepi,qis[i])
        iso_qs = {labels[i]: qis[i] for i in range(len(labels))}
        iso_qs["EEP"] = eepi
        if supplement:
            iso = Isochrone(age=t,
                            eep_name='EEP',
                            mini_name=self.iset.isos[0].mini_name,
                            lteff_name=self.iset.isos[0].lteff_name,
                            llbol_name=self.iset.isos[0].llbol_name,
                            qs=iso_qs)
        else:
            iso = Isochrone(age=t,
                            eep_name='EEP',
                            mini_name=None,
                            lteff_name=None,
                            llbol_name=None,
                            qs=iso_qs)
        return iso

    def _interp_iso_quantity_mass(
        self,
        mini: Quantity["mass"],
        t: Quantity["time"],
        label: str,
        method: str = "pchip",
    ) -> np.float64:
        """
        Uses _construct_iso_isochrone to interpolate a provided quantity in both
        that quantity and in initial mass as a funciton of EEP and then returns the
        value of the quantity at the requested initial mass
        Args:
            mini: the initial mass of the star.
            t: the age of the isochrone.
            label: the label of the quantity to be interpolated.
            method: the interpolation method to use, either pchip or linear
        Returns:
            q_res: the specified quantity at the requested initial mass.
        """
        # construct isochrone for mass/luminosity relationship
        labels = [self.iset.isos[0].mini_name, label]
        iso = self._construct_iso_isochrone(t, labels,method=method,supplement=False)
        qi = iso.qs[label]
        massi = iso.qs[self.iset.isos[0].mini_name]

        mini = mini.to(u.Msun).value
        q_res = pchip_interpolate(massi, qi, mini)
        # make sure values are zero outside the range of masses used
        mask = np.logical_or(mini < min(massi), mini > max(massi))
        q_res = np.ma.masked_array(q_res, mask=mask)
        return q_res

    ######## EEP INTERPOLATION FUNCTIONS ########
    def _get_eep_relation(self, eep:int, labels:list[str]):
        """
        For a given EEP, returns the age, mass and a list of quantities varying
        across all mass tracks at that EEP.
        Args:
            eep: the equivalent evolutionary phase number
            labels: the list of quantity labels to get values for
        Returns:
            age_set: the ages at the given EEP across all tracks
            mass_set: the initial masses at the given EEP across all tracks
            q_set: the list of quantities at the given EEP across all tracks
        """
        nq = len(labels)
        (age_set, mass_set, q_set) = ([], [], [[] for _ in range(nq)])
        for (j,track) in enumerate(self.tset.tracks):
            ages = track.qs[track.age_name]
            if eep <= len(ages):
                # get the age at the given EEP
                age_set.append(ages[eep-1])
                mass_set.append(self.tset.masses[j].value)
                for (i,label) in enumerate(labels):
                    q_set[i].append(track.qs[label][eep-1])
        age_set = np.array(age_set)
        mass_set = np.array(mass_set)*self.tset.masses.unit
        q_set = np.array(q_set)
        return age_set, mass_set, q_set

    def _construct_eep_isochrone(self, t:Quantity["time"], labels:list[str],
                                 method:str="pchip",supplement=True) -> Isochrone:
        """
        Follows the instructions of the MIST0 and MIST1 papers by looping over
        all EEPs and
            1. Looping over all Mass tracks and getting the age at the given EEP
               This then constitutes an age vs. mass relationship which is enforced to be
               monotonic (decreasing as a function of mass)
            2. We interpolate this relationship to find M(t) for the given EEP
            3. We then interpolate the requested quantity, given by "label" in mass for
               the given EEP.
        Args:
            t: The age of the isochrone to be generated, should be a single value.
            labels: The labels of the quantities to be interpolated.
            method: the interpolation method to use, either pchip or linear
            supplement: Whether to supplement the labels with basic quantities
        Returns:
            Isochrone: an isochrone data structure at the requested age t
        """
        # pre-process list of labels to include teff and lbol if wanted
        if supplement:
            labels = self.supplement_labels(labels)

        nq = len(labels)
        interp = self._get_interpolator(method)
        age = array_utils.make_scalar_quantity(t, unit=u.yr).value
        (eeps_,ms_,qs_) = ([],[],[[] for _ in range(nq)])
        for eep in range(1, self.tset.max_eep+1):
            (age_set, mass_set, q_set) = self._get_eep_relation(eep, labels)
            age_test = min(age_set) < age < max(age_set)
            if ((len(age_set) > 0) and age_test):
                eeps_.append(eep)
                age_set = array_utils.make_monotonic_decreasing(mass_set, age_set)
                m = interp(age, age_set[::-1], mass_set[::-1])
                ms_.append(m)
                for i in range(nq):
                    qs_[i].append(interp(m, mass_set, q_set[i]))
        # this is the constructed isochrone for property qs given by "labels"
        iso_qs = {labels[i]: np.array(qs_[i]) for i in range(nq)}
        iso_qs["EEP"] = np.array(eeps_)
        iso_qs["initial_mass"] = np.array(ms_)
        if supplement:
            iso = Isochrone(age=t,
                            eep_name='EEP',
                            mini_name='initial_mass',
                            lteff_name=self.tset.tracks[0].lteff_name,
                            llbol_name=self.tset.tracks[0].llbol_name,
                            qs=iso_qs)
        else:
            iso = Isochrone(age=t,
                            eep_name='EEP',
                            mini_name='initial_mass',
                            lteff_name=None,
                            llbol_name=None,
                            qs=iso_qs)
        return iso

    def _interp_eep_quantity(self, mini:Quantity["mass"], t:Quantity["time"],
                             label:str, method:str="pchip") -> np.float64:
        """
        Uses _construct_eep_isochrone to construct an isochrone at the requested age
        processes the isochrones to insure monotonicity and then interpolates to the
        requested input masses
        Args:
            mini: the initial mass of the star.
            t: the age of the isochrone.
            label: the label of the quantity to be interpolated.
            method: the interpolation method to use, either pchip or linear
        Returns:
            qi: the specified quantity at the requested age (t) and range of initial
                masses (mini)
        """
        interp = self._get_interpolator(method)
        iso = self._construct_eep_isochrone(t, [label], method=method, supplement=False)
        eeps_ = iso.qs[iso.eep_name]
        ms_ = iso.qs[iso.mini_name]
        qs_ = iso.qs[label]
        # make sure masses are monotonically increasing (not guaranteed by above interp)
        ms_ = array_utils.make_monotonic_increasing(eeps_, ms_)
        # finally we interpolate the requested quantity to the requested masses
        mini = mini.to(u.Msun).value
        qi = interp(mini, ms_, qs_)
        # make sure values are zero outside the range of masses used
        mask = np.logical_or(mini < min(ms_), mini > max(ms_))
        qi = np.ma.masked_array(qi, mask=mask)
        return qi

    ######## OTHER HELPER FUNCTIONS ########
    def _get_mmax_age_interp(self):
        """
        Returns the maximum mass as a function of age,
        properly interpreted from either the isochrone or EEP/track data.
        """
        if self.interp_op == "iso":
            lages = self.iset.lages
            mmax = np.array([np.max(iso.qs[iso.mini_name]) for iso in self.iset.isos])
        else:
            # interpolating from EEPs
            lages = np.log10((self.tset.max_ages).to(u.yr).value)[::-1]
            mmax = self.tset.masses[::-1]
            sel = array_utils.index_monotonic(lages)
            lages = lages[sel]
            mmax = mmax[sel]
        return (lages, mmax)

    def _interp_quantity(
        self,
        mini: Quantity["mass"],
        t: Quantity["time"],
        label: str,
        method: str = "pchip",
    ) -> np.float64:
        """
        Helper function to decide between which inerpolation method to use
        and properly format the arguments
        Args:
            mini: the initial mass of the star.
            t: the age of the isochrone. Should be a single time (for now...)
            label: the label of the quantity to be interpolated.
            method: the interpolation method to use, either pchip or linear
        Returns:
            q_res: The quantity specified by label at the requested initial mass.
                   This is masked based on the range of initial masses available
                   at the specified age.
        """

        if np.isscalar(mini.value):
            mini = np.array([mini.value]) * mini.unit

        t = array_utils.make_scalar_quantity(t)
        if self.interp_op == "iso":
            min_age = np.power(10,min(self.iset.lages))*u.yr
            max_age = np.power(10,max(self.iset.lages))*u.yr
        else:
            min_age = min(self.tset.min_ages)
            max_age = max(self.tset.max_ages)
        if (t < min_age) or (t > max_age):
            emsg = f"t = {t.value} yrs is outside the range of applicability: " + \
                    f"{min_age.value} - {max_age.value} yrs"
            raise ValueError(emsg)

        if self.interp_op == "iso":
            # interpolate from isochrones
            q_res = self._interp_iso_quantity_mass(mini, t, label, method=method)
        else:
            # interpolate from EEPs
            q_res = self._interp_eep_quantity(mini, t, label, method=method)

        return q_res

    ######## MAIN PUBLIC OUTPUT FUNCTIONS ########
    def mmax(self, t: Quantity["time"]) -> Quantity["mass"]:
        """
        get the maximum mass of the stellar population that hasn't
        died yet (in e.g. a SN) as a funciton of age, using a cubic spline
        based on maximum mass reported in the isochrone data
        Args:
            t: the time at which to evaluate the derivative, can be an array
        Returns:
            mmax: the initial mass of the most massive star still alive at time t
        """
        if np.isscalar(t.value):
            t = np.array([t.value]) * t.unit

        lt = np.log10(t.to(u.yr).value)

        (lages, mmax) = self._get_mmax_age_interp()
        # cubic spline interpolation in log(age)
        interp = pchip_interpolate(lages,mmax, lt)
        if self.interp_op == "iso":
            max_mass = (self.iset.max_mass).to(u.Msun).value
        else:
            max_mass = np.max(self.tset.masses.to(u.Msun).value)
        interp[np.where(lt < min(lages))] = max_mass
        interp[np.where(lt > max(lages))] = 0.0
        return interp * u.Msun

    def mmaxdot(self, t: Quantity["time"]) -> Quantity["mass"]:
        """
        get the rate of change of the maximum mass of the stellar population
        with respect to time. Uses a cubic spline and takes the derivative
        Args:
            t: the time at which to evaluate the derivative, can be an array
        Returns:
            mmaxdot: the rate at which the maximum mass is changing with respect to time
        """
        if np.isscalar(t.value):
            t = np.array([t.value]) * t.unit

        lt = np.log10(t.to(u.yr).value)

        (lages, mmax) = self._get_mmax_age_interp()
        # return the first derivative of the cubic spline
        unitfac = u.Msun / t.to(u.Myr) / np.log(10)
        interp = pchip_interpolate(lages, mmax, lt, der=1) * unitfac
        interp *= np.logical_and(lt > min(lages), lt < max(lages))
        return interp

    def lbol(
        self, mini: Quantity["mass"], t: Quantity["time"], method: str = "pchip"
    ) -> Quantity["power"]:
        """
        get the bolometric luminosity of a star of initial mass mini at age t
        Args:
            mini: the initial mass of the star. Can be an array
            t: the age of the isochrone. Should be a single time (for now...)
            method: the interpolation method to use, either pchip or linear
        Returns:
            Quantity["power"]: the bolometric luminosity of the star.
        """
        if self.interp_op == "iso":
            llbol_label = self.iset.isos[0].llbol_name
        else:
            llbol_label = self.tset.tracks[0].llbol_name
        logLbol_res = self._interp_quantity(mini, t, llbol_label, method=method)
        return masked_power(10, logLbol_res)*u.Lsun

    def teff(
        self, mini: Quantity["mass"], t: Quantity["time"], method: str = "pchip"
    ) -> Quantity["temperature"]:
        """
        get the atmospheric effective temperature of a star of initial mass mini
        at age t
        Args:
            mini: the initial mass of the star. Can be an array
            t: the age of the isochrone. Should be a single time (for now...)
            method: the interpolation method to use, either pchip or linear
        Returns:
            Quantity["temperature"]: the effective surface temperature of the star.
        """
        if self.interp_op == "iso":
            lteff_label = self.iset.isos[0].lteff_name
        else:
            lteff_label = self.tset.tracks[0].lteff_name
        # interpolating from EEPs
        logTeff_res = self._interp_quantity(mini, t, lteff_label, method=method)
        return masked_power(10, logTeff_res)*u.K
    
    def construct_isochrone(self, t: Quantity["time"],
                            method:str="pchip", labels=None) -> Isochrone:
        """
        Constructs an isochrone at age t with the requested labels
        
        Args:
            t: the age of the isochrone. Should be a single time
            method: the interpolation method to use, either pchip or linear
                    default is pchip (a monotonicity-preserving cubic spline)
            labels: The labels of the quantities to be calcaulated for the isochrone
                    default is an empty list which will return initial mass, log(teff),
                    and log(lbol) only.
        Returns:
            Isochrone: The isochrone object containing initial masses,
                effective temperatures, and bolometric luminosities
                of the stars in the isochrone, along with any other requested quantities.
        """
        t = array_utils.make_scalar_quantity(t)
        if labels is None:
            labels = []
        if self.interp_op == "iso":
            iso = self._construct_iso_isochrone(t, labels, method=method)
        else:
            iso = self._construct_eep_isochrone(t, labels,method=method)
        return iso
