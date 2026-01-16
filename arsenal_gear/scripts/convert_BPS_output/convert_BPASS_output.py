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

def convert_singles(BPASS_directory, output_directory='./arsenal_BPASS', overwrite=False):

    # Create directory if it does not already exists
    Path(output_directory).mkdir(parents=True, exist_ok=True)

    model_directory = BPASS_directory + '/NEWSINMODS/'

    for metal_directory in os.listdir(model_directory):
        
        if not metal_directory.startswith("z"):
            continue

        # Create the zero array
        # Number of files depends on the choice of IMF
        num_files = len(os.listdir(model_directory + '/' + metal_directory))
        # Times for time array
        num_times = 302 # from log(age/yr) = 6 to 9, with 10 times more models than the public models and a 0
        times = np.concatenate((np.zeros(1), np.logspace(6, 9, 301)))
        # Save 6 properties: time, mass, luminosity, temperature, radius, and mass loss rate 
        data = np.zeros((num_files, 6, num_times))
        data[:, 0, :] = np.tile(times, (num_files, 1))

        i = 0
        for model in os.listdir(model_directory + '/' + metal_directory):

            if model.startswith("sneplot"):

                _data = np.genfromtxt(model_directory + '/' + metal_directory 
                                  + '/' + model)

                # Replace NaNs by 0s --> What about files without a companion?
                _data = np.nan_to_num(_data)
                # We want to extract the following values
                _times = _data[:, 1]
                for t in range(len(_times)):

                    if t == 0:
                        _mask = np.where(_times[t] >= times)[0]
                    elif (t == (len(_times) - 1)) and (_times[-1] < times[-1]):
                        _mask = np.arange(302)[np.where(_times[t] >= times)[0][-1]:]
                    else:
                        _mask = np.where((_times[t] >= times) & (_times[t-1] < times))[0]

                    for m in _mask:

                        data[i, 1, m] = _data[t, 5]     # mass in MSun
                        data[i, 2, m] = 10**_data[t, 4] # Lsun
                        if t == (len(_times) - 1):
                            data[i, 2, m] *= 0 # must set to 0 after SN
                        data[i, 3, m] = 10**_data[t, 3] # K
                        data[i, 4, m] = 10**_data[t, 2] # Rsun
                        data[i, 5, m] = (_data[t, 39] * u.Msun / (1.989 * u.s)).to(u.Msun / u.yr).value

            i += 1

        # Sort by mass
        _sort = np.argsort(data[:, 1, 0])
        data_to_save = data[_sort, :, :]

        if (metal_directory + '.npy') not in os.listdir(model_directory) or overwrite:
            np.save(output_directory + '/singles_' + metal_directory + '.npy', data_to_save)
        else:
            print("Cannot save model. Try setting overwrite=True...")

    return

def convert_binaries(BPASS_directory, output_directory='./arsenal_BPASS', overwrite=False):

    # Create directory if it does not already exists
    Path(output_directory).mkdir(parents=True, exist_ok=True)

    model_directory = BPASS_directory + '/NEWBINMODS/NEWBINMODS/'

    for metal_directory in os.listdir(model_directory):
        
        if not metal_directory.startswith("z"):
            continue

        # Create the zero array
        # Number of files depends on the choice of IMF
        num_files = len(os.listdir(model_directory + '/' + metal_directory))
        # Times for time array
        num_times = 502 # from log(age/yr) = 4 to 9, with 10 times more models than the public models and a 0
        times = np.concatenate((np.zeros(1), np.logspace(4, 9, 501)))
        # Save 12 properties: time, primary mass, period, kick (does not evolve), companion mass,
        # luminosities, temperatures, radii, and system mass loss rate
        data = np.zeros((num_files, 12, num_times))
        data[:, 0, :] = np.tile(times, (num_files, 1))

        i = 0
        for model in os.listdir(model_directory + '/' + metal_directory):

            if model.startswith("sneplot"):

                merger = False
                supernova = False

                _data = np.genfromtxt(model_directory + '/' + metal_directory 
                                  + '/' + model)

                # Replace NaNs by 0s --> What about files without a companion?
                _data = np.nan_to_num(_data)
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
                        if (_data[t, 5] > _data[t-1, 5]) and (_data[t, 37] == _data[t-1, 37]) and (merger == False):
                            # Assume merger if M1 increased but M2 remained fixed
                            merger = True

                    if supernova:
                        print("Recomputing orbital elements post-SN...")
                        M1 = _data[t, 5] * u.Msun
                        M2 = _data[t, 37] * u.Msun
                        MR = _data[t, 29] * u.Msun
                        a = 10**_data[t, 35] * u.Rsun
                        v = apply_kick(M1, M2, MR, a, kick_velocity=100 * u.km/u.s, kick_direction=np.array([0, -1, 0]))

                    else:
                        
                        for m in _mask:

                            data[i, 1, m] = _data[t, 5]      # mass in MSun
                            data[i, 2, m] = _data[t, 34]     # period in yr
                            # skip kick for now
                            data[i, 3, m] = _data[t, 29]     # remnant mass in MSun
                            data[i, 4, m] = 10**_data[t, 4]  # Lsun
                            #if t == (len(_times) - 1):
                            #    data[i, 4, m] *= 0 # must set to 0 after SN
                            data[i, 5, m] = 10**_data[t, 3]  # K 
                            data[i, 6, m] = 10**_data[t, 2]  # Rsun
                            # Companion
                            if not merger:
                                data[i, 7, m] = _data[t, 37]     # companion mass in MSun
                                data[i, 8, m] = 10**_data[t, 48] # Lsun, companion
                                data[i, 9, m] = 10**_data[t, 47] # K, companion
                                data[i, 10, m] = 10**_data[t, 46] # Rsun, companion
                            # System mass loss rate
                            data[i, 11, m] = ((_data[t, 39] + _data[t, 40] + _data[t, 43] + _data[t, 44] - _data[t, 41] 
                                               - _data[t, 42]) * u.Msun / (1.989 * u.s)).to(u.Msun / u.yr).value


            i += 1

        # Sort by mass
        _sort = np.argsort(data[:, 1, 0])
        data_to_save = data[_sort, :, :]

        if (metal_directory + '.npy') not in os.listdir(model_directory) or overwrite:
            np.save(output_directory + '/binaries_' + metal_directory + '.npy', data_to_save)
        else:
            print("Cannot save model. Try setting overwrite=True...")

    return
        