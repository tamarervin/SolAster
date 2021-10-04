#!/usr/bin/env python
"""
Tamar Ervin
Date: June 30, 2021

creation of Fig 6 in Milbourne 2019
fraction of observed solar active regions
as function of region area and latitude

"""

import os

import sys

sys.path.append('/Users/tervin/NEID_Solar_analysis')

import time
import datetime
import numpy as np
import pandas as pd

import astropy.units as u

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.coordinates import frames

from skimage.measure import label, regionprops

import tamar.tools.solar_funcs as sfuncs
import tamar.tools.coord_funcs as ctfuncs
import tamar.tools.lbc_funcs as lbfuncs
import tamar.tools.utilities as utils

start_time = time.time()


# name of csv file to store calculations
csv_name = '/Users/tervin/csv_files/milbourne/lat_area_plot.csv'

# print out csv title
print("Beginning calculation of values for csv file: " + csv_name)

# List of header strings
row_contents = ['date_obs', 'co_lat', 'area']

# Append a list as new line to an old csv file
# utils.append_list_as_row(csv_name, row_contents)

# get hmi data products
# cadence = a.Sample(24 * u.hour)  # querying cadence
# start_date = datetime.datetime(2011, 9, 30, 12, 0, 0)
# end_date = datetime.datetime(2011, 12, 7, 12, 0, 0)
time_range = datetime.timedelta(seconds=22)
physobs_list = [a.Physobs.los_velocity, a.Physobs.los_magnetic_field, a.Physobs.intensity]

# create list of dates
# dates_list = [start_date + datetime.timedelta(days=d) for d in range((end_date - start_date).days)]

# create pandas dataframe
data_csv = '/Users/tervin/csv_files/milbourne/milbourne_data.csv'
data_df = pd.read_csv(data_csv)
dates_list = data_df.date_obs.values

for date in dates_list[0:50]:
    # convert the date to a string -- required for use in csv file
    date_str, date_obj = utils.get_dates(date)

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
        save_vals = [date_str, 'not three good files']
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
            save_vals = [date_str, 'missing data product']
            utils.append_list_as_row(csv_name, save_vals)
            pass

        else:
            # transformation for magnetic maps
            mtheta, mphi, mmu = ctfuncs.coord_trans(mmap, R0=mmap.meta['rsun_obs'] * 760e3 / mmap.meta['rsun_ref'], outside_map_val=np.nan)

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

            # calculate magnetic cutoff value
            br_cutoff = Br_cutoff / mmu

            # get active regions
            brthresh = np.full(shape=map_mag_cor.data.shape, fill_value=np.nan)
            use_indices = np.logical_and(mmu > mu_cutoff, np.abs(map_mag_cor.data) > br_cutoff)

            # find isolated pixels
            # convert to integer array for skimage use
            y = use_indices.astype(int)
            # get area
            y_labeled = label(y, connectivity=2, background=0)
            y_area = [props.area for props in regionprops(y_labeled)]

            # get latitude
            mean_lat = np.rint([props.centroid for props in regionprops(y_labeled)]).astype(int)
            theta = np.arccos(mmu)
            lat = [mtheta[l[0], l[1]] for l in mean_lat]

            # get pixel area corresponding to one uHem
            r_sun = mmap.fits_header['rsun_obs']
            pix_dim = mmap.fits_header['cdelt1']
            sun_pix_area = np.pi * (r_sun ** 2) / (pix_dim ** 2)
            hem_area = sun_pix_area / 2 * 10e-6

            # get areas corresponding to latitudes
            arr = [0.09847, -0.09950, 0.06849]
            area_cuts = [(arr[0]*l**2 + arr[1]*l + arr[2]) * hem_area for l in np.abs(lat)]

            # area constraint
            good_area = np.where(np.array(y_area) > area_cuts)
            areas = sorted(good_area[0])
            for use_area in areas:
                val = use_area + 1
                active_indices = np.logical_and(y_labeled == val, y_area[use_area] > area_cuts[use_area])
                lat = np.mean(mtheta[active_indices])
                area = y_area[use_area]
                co_lat = np.pi/2 - lat

                # make array of what we want to save
                save_vals = [co_lat, area]

                # round stuff
                rounded = np.around(save_vals, 3)
                round_vals = [str(mmap.meta['date-obs'])]
                for val in rounded:
                    round_vals.append(val)

                # append these values to the csv file
                utils.append_list_as_row(csv_name, round_vals)

            # print that the date is completed
            print('\nCalculations and save to file complete for ' + date_str)

# print elapsed time
end_time = time.time()

print((end_time - start_time) / 60)


# get latitude
# mean_lat = np.rint([props.centroid for props in regionprops(y_labeled)]).astype(int)
# lat = [theta[l[0], l[1]] for l in mean_lat]
#
# # get pixel area corresponding to one uHem
# r_sun = mmap.fits_header['rsun_obs']
# pix_dim = mmap.fits_header['cdelt1']
# sun_pix_area = np.pi * (r_sun ** 2) / (pix_dim ** 2)
# hem_area = sun_pix_area/2 * 10e-6
#
# # get areas corresponding to latitudes
# area_cuts = [(arr[0] * l ** 2 + arr[1] * l + arr[2])*hem_area for l in np.abs(lat)]
#
#
