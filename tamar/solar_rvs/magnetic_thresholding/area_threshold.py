#!/usr/bin/env python
"""
Tamar Ervin
Date: July 13, 2021

area thresholding of active regions to build
plots comparing active region area to latitude

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

start_time = time.time()

# name of csv file to store calculations
csv_save = '/Users/tervin/csv_files/neid/calcs/area_mag.csv'
csv_name = '/Users/tervin/csv_files/neid/calcs/rvs.csv'

# name of csv to store bad dates
bad_dates_file = '/Users/tervin/csv_files/neid/long_bad_dates.csv'

# print out csv title
print("Beginning calculation of values for csv file: " + csv_name)

# # List of header strings
row_contents = ['date_obs', 'date_jd', 'f_small', 'f_large', 'f_network', 'f_plage', 'f_nonconv',
                'quiet_flux', 'ar_flux', 'conv_flux', 'pol_flux', 'pol_conv_flux', 'vconv_quiet', 'vconv_large',
                'vconv_small']

# Append a list as new line to an old csv file
utils.append_list_as_row(csv_save, row_contents)

# get hmi data products
time_range = datetime.timedelta(seconds=22)
physobs_list = [a.Physobs.los_velocity, a.Physobs.los_magnetic_field, a.Physobs.intensity]

# get dates
df = pd.read_csv(csv_name)
jd_dates = df.date_jd.values

for i, date in enumerate(jd_dates):

    # convert date to usable format
    # date = np.float(date)
    date_str, date_obj, date_jd = utils.get_dates(date)

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
            # get cartesian scales from map
            x, y = ctfuncs.get_map_scales(vmap)

            # get mu and radius array
            mu, r = ctfuncs.get_mu(x, y, vmap)

            # coordinate transformation for maps
            wij, nij, rij = ctfuncs.coord_trans(x, y, r, vmap)

            # remove bad mu values
            vmap, mmap, imap = ctfuncs.fix_mu(mu, [vmap, mmap, imap])

            # calculate relative positions
            deltaw, deltan, deltar, dij = sfuncs.rel_positions(wij, nij, rij, vmap)

            # calculate spacecraft velocity
            vsc = sfuncs.spacecraft_vel(deltaw, deltan, deltar, dij, vmap)

            # optimized solar rotation parameters
            a1 = 14.713
            a2 = -2.396
            a3 = -1.787
            a_parameters = [a1, a2, a3]

            # calculation of solar rotation velocity
            vrot = sfuncs.solar_rot_vel(wij, nij, rij, deltaw, deltan, deltar, dij, vmap, a_parameters)

            # calculate corrected velocity
            corrected_vel = vmap.data - np.real(vsc) - np.real(vrot)

            # corrected velocity maps
            map_vel_cor = sfuncs.corrected_map(corrected_vel, vmap, map_type='Corrected-Dopplergram',
                                               frame=frames.HeliographicCarrington)

            # limb brightening
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

            # radial magnetic data map
            map_mag_cor = sfuncs.corrected_map(Br, mmap, map_type='Corrected-Magnetogram',
                                               frame=frames.HeliographicCarrington)

            ### calculate magnetic threshold
            # magnetic threshold value (G) from Yeo et. al. 2013
            Br_cutoff = 24
            # mu cutoff value
            mu_cutoff = 0.3

            # calculate magnetic threshold
            active, quiet = sfuncs.mag_thresh(mu, mmap, Br_cutoff=Br_cutoff, mu_cutoff=mu_cutoff)

            # calculate intensity threshold
            fac_inds, spot_inds = sfuncs.int_thresh(map_int_cor, active, quiet)

            ### filling factor
            # calculate filling factor
            f_bright, f_spot, f = sfuncs.filling_factor(mu, mmap, active, fac_inds, spot_inds)

            ### calculate the area filling factor
            pixA_hem = ctfuncs.pix_area_hem(wij, nij, rij, vmap)
            area = sfuncs.area_calc(active, pixA_hem)
            f_small, f_large, f_network, f_plage, f_nonconv = sfuncs.area_filling_factor(active, area, mu, mmap, fac_inds)

            ### get the unsigned flux
            quiet_flux, ar_flux, conv_flux, pol_flux, pol_conv_flux = sfuncs.area_unsigned_flux(map_mag_obs, imap, area, active)

            ### get area weighted convective velocities
            vconv_quiet, vconv_large, vconv_small = sfuncs.area_vconv(map_vel_cor, imap, active, area)

            ### add to csv
            # make array of what we want to save
            save_vals = [f_small, f_large, f_network, f_plage, f_nonconv, quiet_flux, ar_flux, conv_flux, pol_flux,
                         pol_conv_flux, vconv_quiet, vconv_large, vconv_small]

            # round stuff
            rounded = np.around(save_vals, 3)
            round_vals = [date_str, date_jd]
            for val in rounded:
                round_vals.append(val)

            # append these values to the csv file
            utils.append_list_as_row(csv_save, round_vals)

    print('Finished calculations for date:', date_str, 'at index:', i)