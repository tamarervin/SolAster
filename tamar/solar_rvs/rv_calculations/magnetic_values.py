#!/usr/bin/env python
"""
Tamar Ervin
Date: June 29, 2021

calculation of magnetic flux and filling factor over long periods
of time

could be interesting for comparisons with AIA images and such
"""

import os

import sys

sys.path.append('/Users/tervin/NEID_Solar_analysis')

import glob
import time
import datetime
import pandas as pd
import numpy as np

import astropy.units as u
from astropy.time import Time
from astropy.io import fits

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.coordinates import frames

import tamar.tools.solar_funcs as sfuncs
import tamar.tools.lbc_funcs as lbfuncs
import tamar.tools.coord_funcs as ctfuncs
import tamar.tools.utilities as utils
from tamar.tools.settings import CsvDir

# start time
start_time = time.time()

# dates
start_date = datetime.datetime(2021, 7, 10, 12, 00, 0)
end_date = datetime.datetime(2021, 8, 6, 12, 00, 0)
dates = [start_date + datetime.timedelta(days=d) for d in range((end_date - start_date).days)]

# name of csv file to store calculations
csv_name = os.path.join(CsvDir.NEID_CALC, 'magnetic_values.csv')

# name of csv to store bad dates
bad_dates_file = os.path.join(CsvDir.NEID_BAD_DATES, 'bad_dates.csv')

# # path to text file
# neid_txt_path = os.path.join(CsvDir.CSV_DIR, 'solar_noon_files_short_name.txt')
#
# file = open(neid_txt_path, 'r')
# lines = file.readlines()
# file.close()

# days_list = [os.path.splitext(line)[0][-15:] for line in lines]
# #
# # # get pieces for date
# year = [d[0:4] for d in days_list]
# month = [d[4:6] for d in days_list]
# day = [d[6:8] for d in days_list]
# hour = [d[9:11] for d in days_list]
# minute = [d[11:13] for d in days_list]
# second = [d[-2:] for d in days_list]
# dates = [d[0:8] for d in days_list]
# dates = np.unique(dates)

# make dates list
# dates_list = [datetime.datetime(*dt) for dt in (map(int, v) for v in zip(year, month, day, hour, minute, second))]

# path to fits files
# fits_path = [os.path.join(CsvDir.NEID_SOLAR, f) for f in lines]

# print out csv title
print("Beginning calculation of values for csv file: " + csv_name)

# List of header strings
# row_contents = ['date_obs', 'date_jd', 'rv_sun', 'v_quiet', 'v_disc', 'v_phot', 'v_conv', 'f_bright',
#                 'f_spot', 'f', 'Bobs']
row_contents = ['date_obs', 'date_jd', 'f_bright', 'f_spot', 'f', 'Bobs']

# Append a list as new line to an old csv file
# utils.append_list_as_row(csv_name, row_contents)

# get hmi data products
time_range = datetime.timedelta(seconds=22)
physobs_list = [a.Physobs.los_velocity, a.Physobs.los_magnetic_field, a.Physobs.intensity]

for i, date in enumerate(dates[-4:]):
    # print(date)
    # file = os.path.join(CsvDir.NEID_SOLAR, date, 'level2', date)
    # spec_fits_files = [i for i in glob.glob(os.path.join(file, '*.fits'))]
    # time = date[0:4] + "-" + date[4:6] + '-' + date[6:8] + 'T12:00:00.00'
    date_str, date_obj, date_jd = utils.get_dates(date)

    # rv_in_file = []
    # for f in spec_fits_files:
    #     main_header = fits.getheader(f)
    #     date_jd = fits.getheader(f)['OBSJD']
    #
    #     if main_header['DRIFTFUN'] == 'simultaneous' and main_header['WAVECAL'] == 'LFCplusThAr':
    #         # get information from fits header
    #         rv_in_file.append(fits.getheader(f, 'CCFS')['CCFRVMOD'])

    # if len(rv_in_file) != 0:
    #     # get the average (binned) rv value
    #     rv_bin = np.average(rv_in_file)

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
        print('Not three good files.')
        save_vals = [date_str, 'not three good files']
        utils.append_list_as_row(bad_dates_file, save_vals)

        pass
    else:
        # convert to map sequence
        map_seq = sunpy.map.Map(sorted(good_files))

        # check for missing data types
        missing_map = False
        # split into data types
        for j, map_obj in enumerate(map_seq):
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
            # append these values to the csv file
            save_vals = [date_str, 'missing data product']
            utils.append_list_as_row(bad_dates_file, save_vals)
            pass

        else:
            # coordinate transformation for maps
            px, py, p, r, d, mu = ctfuncs.coordinates(vmap)
            rw_obs, rn_obs, rr_obs = ctfuncs.vel_coords(px, py, p, r, vmap)

            # remove bad mu values
            vmap, mmap, imap = ctfuncs.fix_mu(mu, [vmap, mmap, imap])

            # # limb brightening
            Lij = lbfuncs.limb_polynomial(imap)

            # calculate corrected data
            Iflat = imap.data / Lij

            # corrected intensity maps
            map_int_cor = sfuncs.corrected_map(Iflat, imap, map_type='Corrected-Intensitygram',
                                               frame=frames.HeliographicCarrington)

            # magnetic noise level
            B_noise = 8

            # calculate unsigned field strength
            Bobs, Br = sfuncs.mag_field(mu, mmap, B_noise)

            # corrected observed magnetic data map
            map_mag_obs = sfuncs.corrected_map(Bobs, mmap, map_type='Corrected-Magnetogram',
                                               frame=frames.HeliographicCarrington)

            # # radial magnetic data map
            # map_mag_cor = sfuncs.corrected_map(Br, mmap, map_type='Corrected-Magnetogram',
            #                                    frame=frames.HeliographicCarrington)

            ### calculate magnetic threshold
            # magnetic threshold value (G) from Yeo et. al. 2013
            Br_cutoff = 24
            # # mu cutoff value
            mu_cutoff = 0.3

            # calculate magnetic threshold
            active, quiet = sfuncs.mag_thresh(mu, mmap, Br_cutoff=Br_cutoff, mu_cutoff=mu_cutoff)

            # calculate intensity threshold
            fac_inds, spot_inds = sfuncs.int_thresh(map_int_cor, active, quiet)

            ### filling factor
            # calculate filling factor
            f_bright, f_spot, f = sfuncs.filling_factor(mu, mmap, active, fac_inds, spot_inds)

            ### unsigned magnetic flux
            # unsigned observed flux
            unsigned_obs_flux = sfuncs.unsigned_flux(map_mag_obs, imap)

            # make array of what we want to save
            # save_vals = [v_quiet, v_disc, v_phot, v_conv, f_bright, f_spot, f, unsigned_obs_flux]
            save_vals = [f_bright, f_spot, f, unsigned_obs_flux]

            # round stuff
            rounded = np.around(save_vals, 3)
            round_vals = [date_str, date_jd]
            for val in rounded:
                round_vals.append(val)

            # append these values to the csv file
            utils.append_list_as_row(csv_name, round_vals)
    # else:
    #     # append these values to the csv file
    #     print('\nNot enough RVs for ' + date_str + ' index: ' + str(i))
    #     utils.append_list_as_row(bad_dates_file, 'Not enough RVs')
    # print that the date is completed
    print('\nCalculations and save to file complete for ' + date_str + ' index: ' + str(i))


# print elapsed time
end_time = time.time()

print((end_time - start_time) / 60)
