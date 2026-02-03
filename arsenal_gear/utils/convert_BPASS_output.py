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

    if ('singles_' + metals + '.h5') not in os.listdir(output_directory) or overwrite:
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

def convert_binaries(BPASS_directory, output_directory, metals='z014', overwrite=False):

    # Create directory if it does not already exists
    Path(output_directory).mkdir(parents=True, exist_ok=True)

    model_directory = 'NEWBINMODS/NEWBINMODS/' + metals

    # Values for disrupted or modified systems, when looking for new systems
    M_single = np.concatenate((np.arange(0.1, 10., 0.1) , 
                               np.arange(10, 100, 1),
                               np.arange(100, 325, 25)))

    # Times for time array
    num_times = 502 # from log(age/yr) = 6 to 9, with 10 times more models than the public models and a 0
    times = np.concatenate((np.zeros(1), np.logspace(4, 9, 501)))

    data = None

    tar = tarfile.open(BPASS_directory + "/bpass-v2.2-newmodels.tar.gz", "r:gz")

    for model_file in tar:

        if model_file.name.startswith(model_directory + '/sneplot'):
            print(model_file.name)

            output_data = np.zeros((1, 6, num_times))

            f = tar.extractfile(model_file)
            input_data = np.genfromtxt(f)

            merger = False
            supernova = False
            rejuvenated = False
            accreted_mass = 0
            accretion_time = None
            effective_time = None
            rejuvenated_file = None

            # Replace NaNs by 0s --> What about files without a companion?
            input_data = np.nan_to_num(input_data)

            # Effective companion mass
            M2_eff = input_data[0, 37]
                
            # We want to extract the following values
            input_times = input_data[:, 1]
            
            for t in range(len(input_times)):

                if t == 0:
                    _mask = np.where(input_times[t] >= times)[0]
                elif (t == (len(input_times) - 1)) and (input_times[-1] < times[-1]):
                    _mask = np.arange(len(times))[np.where(input_times[t] >= times)[0][-1]:]
                    # If initial mass above 8 Msun, assume SN - to adjust? 
                    if input_data[0, 5] > 8:
                        supernova = True
                else:
                    _mask = np.where((input_times[t] >= times) & (input_times[t-1] < times))[0]
                    if (input_data[t, 5] > input_data[t-1, 5]) and (input_data[t, 37] == input_data[t-1, 37]):# and (merger == False):
                        # Assume merger if M1 increased but M2 remained fixed
                        # Note that this can take place over several timesteps
                        merger = True
                        accreted_mass += input_data[t, 5] - input_data[t-1, 5]
                    elif (input_data[t, 37] > input_data[t-1, 37]) and (merger == False):
                        # Assume accretion if M2 increased
                        M2_eff = input_data[t, 37]
                        if (input_data[t, 37]/input_data[0, 37] >= 1.05) and (input_data[0, 37] >= 2):
                            rejuvenated = True
                        if accretion_time == None:
                            accretion_time = input_times[t]

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
                        effective_time = input_times[t] - accretion_time
                    else:
                        effective_time = input_times[t]

                    rejuvenated_file = 'NEWSINMODS/' + metals + '/sneplot-' + metals + '-' + M2_closest


                else:
                    
                    for m in _mask:

                        output_data[0, 0, m] = input_data[t, 5]      # mass in MSun
                        output_data[0, 1, m] = input_data[t, 4]  # Lsun
                        if t == (len(input_times) - 1):
                            data[0, 1, m] *= 0 # must set to 0 after SN
                        output_data[0, 2, m] = input_data[t, 3]  # K 
                        # Companion
                        if not merger:
                            output_data[0, 3, m] = input_data[t, 37]     # companion mass in MSun
                            output_data[0, 4, m] = input_data[t, 48] # Lsun, companion
                            output_data[0, 5, m] = input_data[t, 47] # K, companion
                        if merger: 
                            if (input_data[t, 5] - input_data[t-1, 5]) > 0: # If still accreting
                                output_data[0, 3, m] = input_data[t, 37] - accreted_mass
                                output_data[0, 4, m] = input_data[t, 48] # Lsun, companion
                                output_data[0, 5, m] = input_data[t, 47] # K, companion

                if rejuvenated_file:

                    f = tar.extractfile(rejuvenated_file)
                    M2_data = np.genfromtxt(f)

                    # Replace NaNs by 0s --> What about files without a companion?
                    M2_data = np.nan_to_num(M2_data)

                    M2_times = M2_data[:, 1][M2_data[:, 1] > effective_time]
                        
                    for t in range(len(M2_times)):

                        if (t == (len(M2_times) - 1)) and (M2_times[-1] < times[-1]):
                            _mask = np.arange(len(times))[np.where(M2_times[t] >= times)[0][-1]:]
                        else:
                            _mask = np.where((M2_times[t] >= times) & (M2_times[t-1] < times))[0]


                        for m in _mask:

                            # Keep primary mass
                            output_data[0, 0, m] = output_data[0, 0, m-1]      # mass in MSun
                            # Companion properties
                            output_data[0, 3, m] = M2_data[t, 5]      # mass in MSun
                            output_data[0, 4, m] = M2_data[t, 4]  # Lsun
                            if t == (len(M2_times) - 1):
                                output_data[0, 4, m] *= 0 # must set to 0 after SN
                            output_data[0, 5, m] = M2_data[t, 3]  # K 

            if data is None:
                data = output_data
            else:
                data = np.concatenate((data, output_data))

    tar.close()
        
    # Sort by mass
    _sort = np.argsort(data[:, 0, 0])[::-1]
    data_to_save = data[_sort, :, :]

    if ('binaries_' + metals + '.h5') not in os.listdir(output_directory) or overwrite:
        print("Saving processed data to", output_directory)

        # Times
        times = np.round(times, 2).astype('str')
        for t in range(len(times)):
            times[t] = times[t].ljust(4, '0')

        # Masses as strings
        masses = data_to_save[:, 0, 0].astype('str')
        for m in range(len(masses)):
            masses[m] = masses[m].ljust(4, '0')

        ds = xr.DataArray(data, coords=[("Model", masses), ("Property", ["Mass 1 (MSun)", "log L_bol 1 (LSun)", "log T_eff 1 (K)",
                                                                        "Mass 2 (MSun)", "log L_bol 2 (LSun)", "log T_eff 2 (K)"]), ("Time (log t/yr)", times)])
        ds.to_netcdf(output_directory + '/binaries_' + metals + '.h5')

    else:
        print("Cannot save model. Try setting overwrite=True...")

    return
        