"""
data_converter.py
================================

This file defines the interface to various binary evolution models
through downloading, reorganzing and interpreting their outputs.
"""

import os
from pathlib import Path
from abc import ABC, abstractmethod

import numpy as np
import xarray as xr

from .be_data_structures import SingleStarTable, BinaryStarTable


class BinaryEvolutionConverter(ABC):
    """
    Abstract base class for converting binary evolution models from various
    sources to a standard Arsenal binary evolution format.
    """
    def __init__(self, **kwargs) -> None:
        # [Fe/H]
        self.met = kwargs.get('met', 0.014)
        # Directories to read and write data
        self.input_dir  = kwargs.get('input_dir', None)
        self.output_dir = kwargs.get('output_dir', None)
        # Output times
        self.output_times = kwargs.get('output_times', np.logspace(4, 9, 501))
        self.overwrite = kwargs.get('overwrite', False)

    @abstractmethod
    def convert_single_data(self) -> SingleStarTable:
        """
        Abstract method for converting single star evolutionary tracks
        into an Arsenal stellar evolution table.
        """

    @abstractmethod
    def convert_binary_data(self) -> BinaryStarTable:
        """
        Abstract method for converting binary evolutionary tracks into
        an Arsenal binary evolution table.
        """

class BPASSConverter(BinaryEvolutionConverter):
    """
    Class for reading BPASS data and converting it to the Arsenal binary evolution
    format. This is an instantiation of the BinaryEvolutionConverter base class.
    """
    mets = ["zem5", "zem4", "z001", "z002", "z003", "z004", "z006",
            "z008", "z010", "z014", "z020", "z030", "z040"]

    def __init__(self, **kwargs) -> None:
        """
        Args:
            kwargs: Keyword arguments for the binary evolution table.

        Methods:
            convert_single_data     Processes single stellar track data into a SingleStarTable
            convert_binary_data     Processes binary stellar track data into a BinaryStarTable
        """
        # set input parameters
        super().__init__(**kwargs)

        if self.met>=1e-3:
            self.metstr = "z" + str(int(self.met * 1000)).zfill(3)
        else:
            self.metstr = "zem" + str(-1*int(np.log10(self.met)))
        if self.metstr not in self.mets:
            raise ValueError("Metallicity must be one of: " + str(self.mets))

        # Consistent format for directories
        if self.input_dir[-1] == "/":
            self.input_dir: str = self.input_dir
        else:
            self.intput_dir: str = self.input_dir + "/"
        if self.output_dir[-1] == "/":
            self.output_dir: str = self.output_dir
        else:
            self.output_dir: str = self.output_dir + "/"
        

    def convert_single_data(self):
        """
        Converts BPASS data for single stars into an Arsenal-readable 
        SingleStarTable. 
        TODO (06.02.2026): Make the number of time outputs and the 
        returned properties user parameters.
        """

        # Create directory if it does not already exists
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

        model_directory = self.input_dir + '/singles/' + self.metstr

        # Create the zero array
        # Number of files depends on the choice of mass limits for the IMF
        num_files = len(os.listdir(model_directory)) 
        # Times for time array
        num_times = len(self.output_times) + 1
        times = np.concatenate((np.zeros(1), self.output_times))
        # Save properties as a function of time: (1) mass, (2) bolometric luminosity, 
        # (3) surface temperature, (4) radius
        data = np.zeros((num_files, 4, num_times))

        i = 0
        for model in os.listdir(model_directory):

            if model.startswith("sneplot"):

                _data = np.genfromtxt(model_directory + '/' + model)

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
                        
                        data[i, 0, m] = _data[t, 5] # mass in Msun
                        data[i, 1, m] = _data[t, 4] # log Lbol in Lsun
                        data[i, 2, m] = _data[t, 3] # log Teff in K
                        data[i, 3, m] = _data[t, 2] # log R in Rsun
                        if t == (len(_times) - 1):
                            data[i, 1, m] *= 0      # must set to 0 after SN
                            data[i, 2, m] *= 0

                i += 1

        # Sort by mass
        _sort = np.argsort(data[:, 0, 0])
        data = data[_sort, :, :]


        if ('singles_' + self.metstr + '.h5') not in os.listdir(self.output_dir) or self.overwrite:
            print("Saving processed data to", self.output_dir)

            masses = data[:, 0, 0]

            ds = xr.Dataset(data_vars=dict(mass=(["model", "time"], data[:, 0, :]),
                                           log_Lbol=(["model", "time"], data[:, 1, :]),
                                           log_Teff=(["model", "time"], data[:, 2, :]),
                                           log_R=(["model", "time"], data[:, 3, :])),
                            coords=dict(model=("model", masses), time=("time", times)),
                            attrs=dict(description="BPASS evolution data for single stars at Z=" + str(self.met)))
            ds.to_netcdf(self.output_dir + '/singles_' + self.metstr + '.h5')

        else:
            print("Cannot save model. Try setting overwrite=True...")

        return


    def convert_binary_data(self):
        """
        Converts BPASS data for binary stars into an Arsenal-readable 
        BinaryStarTable. 
        TODO (06.02.2026): Make the number of time outputs and the 
        returned properties user parameters.
        """

        # Create directory if it does not already exists
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

        model_directory = self.input_dir + '/binaries/' + self.metstr

        # Values for disrupted or modified systems, when looking for new systems
        M_single = np.concatenate((np.arange(0.1, 10., 0.1) , 
                                   np.arange(10, 100, 1),
                                   np.arange(100, 325, 25)))

        # Create the zero array
        # Number of files depends on the choice of IMF
        num_files = len(os.listdir(model_directory))
        # Times for time array
        num_times = len(self.output_times) + 1
        times = np.concatenate((np.zeros(1), self.output_times))
        # Save properties as a function of time: (1)  primary mass, 
        # (2) primary bolometric luminosity, (3) primary surface temperature, 
        # (4) primary radius, (5)-(8) for the companion
        data = np.zeros((num_files, 7, num_times))

        M1_for_models = np.zeros(num_files)
        M2_for_models = np.zeros(num_files)
        P_for_models = np.zeros(num_files)

        i = 0
        for model in os.listdir(model_directory):

            if model.startswith("sneplot"):

                merger = False
                supernova = False
                rejuvenated = False
                accreted_mass = 0
                accretion_time = None
                effective_time = None
                rejuvenated_file = None

                _data = np.genfromtxt(model_directory + '/' + model)

                # Replace NaNs by 0s --> What about files without a companion?
                _data = np.nan_to_num(_data)

                # Effective companion mass
                M2_eff = _data[0, 37]
                M2_for_models[i] = M2_eff

                # Initial primary mass & period
                M1_for_models[i] = _data[0, 5]
                P_for_models[i] = _data[0, 34]
                
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

                        rejuvenated_file = self.input_dir + '/singles/' + self.metstr + '/sneplot-' + self.metstr + '-' + M2_closest


                    else:
                        
                        for m in _mask:

                            data[i, 0, m] = _data[t, 5]  # mass in MSun
                            data[i, 1, m] = _data[t, 4]  # Lsun
                            if t == (len(_times) - 1):
                                data[i, 1, m] *= 0 # must set to 0 after SN
                            data[i, 2, m] = _data[t, 3]  # K 
                            # Companion
                            if not merger:
                                data[i, 3, m] = _data[t, 37] # companion mass in MSun
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
        

        if ('binaries_' + self.metstr + '.h5') not in os.listdir(self.output_dir) or self.overwrite:
            print("Saving processed data to", self.output_dir)


            # Masses as strings
            masses = M1_for_models
            model_masses = np.sort(np.unique(masses))[::-1]

            mass_ratios = np.round(M2_for_models/M1_for_models, 1)
            model_mratios = np.sort(np.unique(mass_ratios))
            model_mratios = model_mratios[model_mratios > 0]

            periods = np.round(np.log10(P_for_models*365.25), 1)
            model_periods = np.sort(np.unique(periods))
            model_periods = model_periods[model_periods < 5]

            model_orbits = (np.outer(model_mratios*1000, np.ones(len(model_periods))) + np.outer(np.ones(len(model_mratios)), model_periods*10)).flatten()
            
            def prop_arr(times, masses, mass_ratios, periods, model_m, model_o, data):

                assert len(masses) == len(mass_ratios)
                assert len(masses) == len(periods)
                assert len(data) == len(masses)

                arr = np.zeros((len(times), len(model_m), len(model_o)))
                for i in range(len(masses)):
                    if (mass_ratios[i] > 0) and (periods[i] < 5):
                        a = np.where(model_m == masses[i])[0][0]
                        b = np.where(model_o == (mass_ratios[i]*1000 + periods[i]*10))[0][0]
                        arr[:, a, b] = data[i]
                return arr

            mass1_arr = prop_arr(times, masses, mass_ratios, periods, 
                                 model_masses, model_orbits, data[:, 0, :])
            mass2_arr = prop_arr(times, masses, mass_ratios, periods, 
                                 model_masses, model_orbits, data[:, 3, :])
            log_Lbol1_arr = prop_arr(times, masses, mass_ratios, periods, 
                                     model_masses, model_orbits, data[:, 1, :])
            log_Lbol2_arr = prop_arr(times, masses, mass_ratios, periods, 
                                     model_masses, model_orbits, data[:, 4, :])
            log_Teff1_arr = prop_arr(times, masses, mass_ratios, periods, 
                                     model_masses, model_orbits, data[:, 2, :])
            log_Teff2_arr = prop_arr(times, masses, mass_ratios, periods, 
                                     model_masses, model_orbits, data[:, 5, :])

            print(mass2_arr)

            ds = xr.Dataset(data_vars=dict(mass1=(["time", "model_m", "model_o"], mass1_arr),
                                           log_Lbol1=(["time", "model_m", "model_o"], log_Lbol1_arr),
                                           log_Teff1=(["time", "model_m", "model_o"], log_Teff1_arr),
                                           mass2=(["time", "model_m", "model_o"], mass2_arr),
                                           log_Lbol2=(["time", "model_m", "model_o"], log_Lbol2_arr),
                                           log_Teff2=(["time", "model_m", "model_o"], log_Teff2_arr),),
                            coords=dict(model_m=(["model_m"], model_masses), model_o=(["model_o"], model_orbits),
                                        time=("time", times)),
                            attrs=dict(description="BPASS evolution data for binary stars at Z=" + str(self.met)))

            ds.to_netcdf(self.output_dir + '/binaries_' + self.metstr + '.h5')

        else:
            print("Cannot save model. Try setting overwrite=True...")

        return
