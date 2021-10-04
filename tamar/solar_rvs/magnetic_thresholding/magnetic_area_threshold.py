#!/usr/bin/env python
"""
Tamar Ervin
Date: June 29, 2021

messing with area constraint of magnetic threshold

"""

import os

import sys

sys.path.append('//')

import time
import datetime
import pandas as pd
import numpy as np

import astropy.units as u
from astropy.time import Time

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.coordinates import frames

import tamar.tools.solar_funcs as sfuncs
import tamar.tools.lbc_funcs as lbc
import tamar.tools.coord_funcs as ctfuncs
import tamar.tools.utilities as utils

start_time = time.time()


# name of csv file to store calculations
csv_name = '/Users/tervin/csv_files/milbourne_mag_area.csv'

# print out csv title
print("Beginning calculation of values for csv file: " + csv_name)

# List of header strings
row_contents = ['date_obs', 'v_quiet', 'v_disc', 'v_phot', 'v_conv', 'f', 'Bobs', 'mag_area']

# Append a list as new line to an old csv file
utils.append_list_as_row(csv_name, row_contents)

# get hmi data products
cadence = a.Sample(24 * u.hour)  # querying cadence
start_date = datetime.datetime(2011, 9, 30, 12, 0, 0)
end_date = datetime.datetime(2011, 12, 7, 12, 0, 0)
time_range = datetime.timedelta(seconds=22)
physobs_list = [a.Physobs.los_velocity, a.Physobs.los_magnetic_field, a.Physobs.intensity]

# create list of dates
# dates_list = [start_date + datetime.timedelta(days=d) for d in range((end_date - start_date).days)]

# read in Haywood 2016 data file (use the csv I made)
# create pandas dataframe
data_csv = '/Users/tervin/csv_files/milbourne_data.csv'
data_df = pd.read_csv(data_csv)

dates_list = data_df.date_obs.values

for date in dates_list[:10]:
    # convert the date to a string -- required for use in csv file
    if isinstance(date, str):
        date_str = date
        date_obj = datetime.datetime.strptime(date, '%Y-%m-%dT%H:%M:%S.%f')
    elif isinstance(date, float):
        t = Time(date, format='jd')
        date_obj = t.datetime
        date_str = date_obj.strftime('%Y-%m-%dT%H:%M:%S')
    else:
        date_obj = date
        date_str = date.strftime('%Y-%m-%dT%H:%M:%S')

    # pull image within specified time range
    result = Fido.search(a.Time(str(date_obj - time_range), str(date_obj + time_range)),
                         a.Instrument.hmi, physobs_list[0] | physobs_list[1] | physobs_list[2])

    # add file to list
    file_download = Fido.fetch(result)

    # remove unusable file types
    good_files = []
    for file in file_download:
        name, extension = os.path.splitext(file)
        if extension == '.fits':
            good_files.append(file)

    if len(good_files) != 3:
        # add the data
        # append these values to the csv file
        save_vals = [str(mmap.meta['date-obs']), 'not three good files']
        utils.append_list_as_row(csv_name, save_vals)

        pass
    else:
        # convert to map sequence
        map_seq = sunpy.map.Map(sorted(good_files))

        # check for missing data types
        missing_map = False
        # split into data types
        for i, map_obj in enumerate(map_seq):
            if map_obj.meta['content'] == 'DOPPLERGRAM':
                vmap = map_obj
            elif map_obj.meta['content'] == 'MAGNETOGRAM':
                mmap = map_obj
            elif map_obj.meta['content'] == 'CONTINUUM INTENSITY':
                imap = map_obj
            else:
                missing_map = True

        if missing_map:
            print("Missing a data product for " + date_str)

            # add the data
            # append these values to the csv file
            save_vals = [str(mmap.meta['date-obs']), 'missing data product']
            utils.append_list_as_row(csv_name, save_vals)
            pass

        else:
            # transformation for velocity maps
            vtheta, vphi, vmu = ctfuncs.coord_trans(vmap, R0=vmap.meta['rsun_obs'] * 760e3 / vmap.meta['rsun_ref'], outside_map_val=np.nan)

            # transformation for intensity maps
            itheta, iphi, imu = ctfuncs.coord_trans(imap, R0=imap.meta['rsun_obs'] * 760e3 / imap.meta['rsun_ref'], outside_map_val=np.nan)

            # transformation for magnetic maps
            mtheta, mphi, mmu = ctfuncs.coord_trans(mmap, R0=mmap.meta['rsun_obs'] * 760e3 / mmap.meta['rsun_ref'], outside_map_val=np.nan)

            # calculate relative positions
            deltaw, deltan, deltar, dij = sfuncs.rel_positions(vtheta, vphi, vmap)

            # calculate spacecraft velocity
            vsc = sfuncs.spacecraft_vel(deltaw, deltan, deltar, dij, vmap)

            # optimized solar rotation parameters
            a1 = 14.713
            a2 = -2.396
            a3 = -1.787
            a_parameters = [a1, a2, a3]

            # calculation of solar rotation velocity
            vrot = sfuncs.solar_rot_vel(a_parameters, vtheta, vphi, deltaw, deltan, deltar, dij, vmap)

            # calculate corrected velocity
            corrected_vel = vmap.data - vrot - vsc

            # corrected velocity maps
            map_vel_cor = sfuncs.corrected_map(corrected_vel, vmap, map_type='Corrected-Dopplergram',
                                               frame=frames.HeliographicCarrington)

            # limb brightening
            Lij = lbc.limb_polynomial(imap)

            # calculate corrected data
            Iflat = imap.data / Lij

            # corrected intensity maps
            map_int_cor = sfuncs.corrected_map(Iflat, imap, map_type='Corrected-Intensitygram',
                                               frame=frames.HeliographicCarrington)

            # magnetic noise level
            B_noise = 8

            # calculate unsigned field strength
            Bobs, Br = sfuncs.mag_field(mmu, mmap, B_noise)

            # corrected observed magnetic data map
            map_mag_obs = sfuncs.corrected_map(Bobs, mmap, map_type='Corrected-Magnetogram',
                                               frame=frames.HeliographicCarrington)

            # radial magnetic data map
            map_mag_cor = sfuncs.corrected_map(Br, mmap, map_type='Corrected-Magnetogram',
                                               frame=frames.HeliographicCarrington)

            ### calculate magnetic threshold
            # magnetic threshold value (G) from Yeo et. al. 2013
            Br_cutoff = 24
            # mu cutoff value
            mu_cutoff = 0.3

            # area constraints list
            mag_area = [20]

            for area in mag_area:
                # calculate magnetic threshold
                Brthresh, mag_weights = sfuncs.mag_thresh(mmu, map_mag_cor, Br_cutoff=Br_cutoff,
                                                                            mu_cutoff=mu_cutoff, area=area)

                # create array of active weights (flipped from mag weights)
                act_weights = sfuncs.active_weights(mag_weights)

                # thresholded magnetic maps
                map_mag_thresh = sfuncs.corrected_map(Brthresh, mmap, map_type='Magnetic-Threshold',
                                                      frame=frames.HeliographicCarrington)

                # calculate intensity threshold
                faculae, sunspots, fac_ind, spot_ind, int_cutoff = sfuncs.int_thresh(map_int_cor, map_mag_cor, mag_weights)

                # create threshold array
                thresh_arr = sfuncs.thresh_map(fac_ind, spot_ind)

                # full threshold maps
                map_full_thresh = sfuncs.corrected_map(thresh_arr, mmap, map_type='Threshold', frame=frames.HeliographicCarrington)

                ### velocity contribution due to convective motion of quiet-Sun
                v_quiet = sfuncs.v_quiet(map_vel_cor, imap, mag_weights)

                ### velocity contribution due to rotational Doppler imbalance of active regions (faculae/sunspots)
                # calculate photospheric velocity
                v_phot = sfuncs.v_phot(mag_weights, act_weights, Lij, vsc, imap, vmap)

                ### velocity contribution due to suppression of convective blueshift by active regions
                # calculate disc-averaged velocity
                v_disc = sfuncs.v_disc(map_vel_cor, imap)

                # calculate convective velocity
                v_conv = v_disc - v_quiet - v_phot

                ### filling factor
                # calculate filling factor
                filling_factor = sfuncs.filling_factor(map_mag_obs, act_weights)

                ### unsigned magnetic flux
                # unsigned observed flux
                unsigned_obs_flux = sfuncs.unsigned_flux(map_mag_obs, imap)

                # make array of what we want to save
                save_vals = [str(mmap.meta['date-obs']), v_quiet,  v_disc, v_phot, v_conv, filling_factor,
                             unsigned_obs_flux, area]

                # append these values to the csv file
                utils.append_list_as_row(csv_name, save_vals)

            # print that the date is completed
            print('Calculations and save to file complete for ' + date_str)

# print elapsed time
end_time = time.time()

print((end_time - start_time) / 60)