"""
This script converts the downloaded BPASS output into an
Arsenal-readable output. We may decide to make either the script
or its output part of the Arsenal distribution.
"""

import os
import astropy.units as u
import numpy as np
from astropy.units import Quantity
from pathlib import Path
import tarfile
import xarray as xr

def velocity_from_orbit(M1, M2, semi_major_axis):
    v = (c.G * (M1+M2) / semi_major_axis)**(1./2)
    return v

def apply_kick(M1, M2, MR, semi_major_axis, kick_velocity, kick_direction=np.array([1, 0, 0])):
    # https://iopscience.iop.org/article/10.1086/340494
    # https://arxiv.org/pdf/2504.16161v2
    v = velocity_from_orbit(M1, M2, semi_major_axis)
    v1 = np.array([0, M2/(M1+M2), 0]) * v
    v2 = np.array([0, -M1/(M1+M2), 0]) * v
    vR = v1 + kick_direction*kick_velocity

    # New angular momentum vector
    ang_mom = np.cross(np.array([1, 0, 0]) * semi_major_axis, vR)
    ecc_vec = (np.cross(vR-v2, ang_mom)/(c.G * (MR + M2))) - np.array([1, 0, 0])

    new_a = np.dot(ang_mom, ang_mom)/(c.G * (MR + M2) * (1 - np.dot(ecc_vec, ecc_vec)))

    new_P = 2*np.pi * (new_a**3/(c.G * (MR + M2)))**(1./2)
    
    return vR

def convert_singles(BPASS_directory, output_directory, metals='z014', overwrite=False):

    # Create directory if it does not already exists
    Path(output_directory).mkdir(parents=True, exist_ok=True)

    model_directory = 'NEWSINMODS/' + metals

    # Times for time array
    num_times = 502 # from log(age/yr) = 6 to 9, with 10 times more models than the public models and a 0
    times = np.concatenate((np.zeros(1), np.logspace(4, 9, 501)))

    data = None

    tar = tarfile.open(BPASS_directory + "/bpass-v2.2-newmodels.tar.gz", "r:gz")

    for model_file in tar:

        if model_file.name.startswith(model_directory + '/sneplot'):
            print(model_file.name)

            output_data = np.zeros((1, 3, num_times))

            f = tar.extractfile(model_file)
            input_data = np.genfromtxt(f)

            # Replace NaNs by 0s --> What about files without a companion?
            input_data = np.nan_to_num(input_data)
            # We want to extract the following values
            input_times = input_data[:, 1]
            for t in range(len(input_times)):

                if t == 0:
                    _mask = np.where(input_times[t] >= times)[0]
                elif (t == (len(input_times) - 1)) and (input_times[-1] < times[-1]):
                    _mask = np.arange(len(times))[np.where(input_times[t] >= times)[0][-1]:]
                else:
                    _mask = np.where((input_times[t] >= times) & (input_times[t-1] < times))[0]

                for m in _mask:

                    output_data[0, 0, m] = input_data[t, 5] # mass in MSun
                    output_data[0, 1, m] = input_data[t, 4] # Lsun
                    if t == (len(input_times) - 1):
                            output_data[0, 1, m] *= 0 # must set to 0 after SN
                    output_data[0, 2, m] = input_data[t, 3] # K

            if data is None:
                data = output_data
            else:
                data = np.concatenate((data, output_data))
        
    tar.close()
        
    # Sort by mass
    _sort = np.argsort(data[:, 0, 0])[::-1]
    data_to_save = data[_sort, :, :]

    if ('singles_' + metals + '.npy') not in os.listdir(output_directory) or overwrite:
        print("Saving processed data to", output_directory)

        # Times
        times = np.round(times, 2).astype('str')
        for t in range(len(times)):
            times[t] = times[t].ljust(4, '0')

        # Masses as strings
        masses = data_to_save[:, 0, 0].astype('str')
        for m in range(len(masses)):
            masses[m] = masses[m].ljust(4, '0')

        ds = xr.DataArray(data, coords=[("Model", masses), ("Property", ["Mass (MSun)", "log L_bol (LSun)", "log T_eff (K)"]), ("Time (log t/yr)", times)])
        ds.to_netcdf(output_directory + '/singles_' + metals + '.h5')

    else:
        print("Cannot save model. Try setting overwrite=True...")

    return

def convert_binaries(BPASS_directory, output_directory='./arsenal_BPASS', overwrite=False):

    # Create directory if it does not already exists
    Path(output_directory).mkdir(parents=True, exist_ok=True)

    model_directory = BPASS_directory + '/NEWBINMODS/NEWBINMODS/'

    # Values for disrupted or modified systems, when looking for new systems
    M_single = np.concatenate((np.arange(0.1, 10., 0.1) , 
                               np.arange(10, 100, 1),
                               np.arange(100, 325, 25)))

    for metal_directory in os.listdir(model_directory):
        
        if not metal_directory.startswith("z"):
            continue

        # Create the zero array
        # Number of files depends on the choice of IMF
        num_files = len(os.listdir(model_directory + '/' + metal_directory)) - 1 # -1 for directory itself
        # Times for time array
        num_times = 502 # from log(age/yr) = 4 to 9, with 10 times more models than the public models and a 0
        times = np.concatenate((np.zeros(1), np.logspace(4, 9, 501)))
        # Save 8 properties: masses, luminosities, temperatures, radii,
        data = np.zeros((num_files, 8, num_times))
        #data[:, 0, :] = np.tile(times, (num_files, 1))

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

                _data = np.genfromtxt(model_directory + '/' + metal_directory 
                                  + '/' + model)

                # Replace NaNs by 0s --> What about files without a companion?
                _data = np.nan_to_num(_data)

                # Effective companion mass
                M2_eff = _data[0, 37]
                
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
                        if accretion_time:
                            effective_time = _times[t] - accretion_time
                        else:
                            effective_time = _times[t]

                        rejuvenated_file = BPASS_directory + '/NEWSINMODS/' + metal_directory + '/sneplot-' + metal_directory + '-' + M2_closest


                    else:
                        
                        for m in _mask:

                            data[i, 0, m] = _data[t, 5]      # mass in MSun
                            data[i, 1, m] = 10**_data[t, 4]  # Lsun
                            if t == (len(_times) - 1):
                                data[i, 1, m] *= 0 # must set to 0 after SN
                            data[i, 2, m] = 10**_data[t, 3]  # K 
                            data[i, 3, m] = 10**_data[t, 2]  # Rsun
                            # Companion
                            if not merger:
                                data[i, 4, m] = _data[t, 37]     # companion mass in MSun
                                data[i, 5, m] = 10**_data[t, 48] # Lsun, companion
                                data[i, 6, m] = 10**_data[t, 47] # K, companion
                                data[i, 7, m] = 10**_data[t, 46] # Rsun, companion
                            if merger: 
                                if (_data[t, 5] - _data[t-1, 5]) > 0: # If still accreting
                                    data[i, 4, m] = _data[t, 37] - accreted_mass
                                    data[i, 5, m] = 10**_data[t, 48] # Lsun, companion
                                    data[i, 6, m] = 10**_data[t, 47] # K, companion
                                    data[i, 7, m] = 10**_data[t, 46] # Rsun, companion

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
                                data[i, 0, m] = data[i, 0, m-1]      # mass in MSun
                                # Companion properties
                                data[i, 4, m] = _data_M2[t, 5]      # mass in MSun
                                data[i, 5, m] = 10**_data_M2[t, 4]  # Lsun
                                if t == (len(_times_M2) - 1):
                                    data[i, 5, m] *= 0 # must set to 0 after SN
                                data[i, 6, m] = 10**_data_M2[t, 3]  # K 
                                data[i, 7, m] = 10**_data_M2[t, 2]  # Rsun


                i += 1


        # Sort by mass
        _sort = np.argsort(data[:, 0, 0])[::-1]
        data_to_save = data[_sort, :, :]

        if (metal_directory + '.npy') not in os.listdir(model_directory) or overwrite:
            np.save(output_directory + '/binaries_' + metal_directory + '.npy', data_to_save)
        else:
            print("Cannot save model. Try setting overwrite=True...")

    return
        