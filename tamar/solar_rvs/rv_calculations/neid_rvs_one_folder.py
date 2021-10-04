#!/usr/bin/env python
"""
Tamar Ervin
Date: August 12, 2021

outline to calculate RVs for NEID dates
- using the one hour averaged around Solar noon
to calculate RV components

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

start_time = time.time()

# name of csv file to store calculations
csv_name = os.path.join(CsvDir.NEID_CALC, 'rvs.csv')

# name of csv to store bad dates
bad_dates_file = os.path.join(CsvDir.NEID_BAD_DATES, 'long_bad_dates.csv')

# get dates from neid files
spec_fits_files = [i for i in glob.glob(os.path.join(CsvDir.NEID_HOUR, '*.fits'))]
dates = [s.split('_')[1] for s in spec_fits_files]
dates = [s.split('.')[0] for s in dates]
days = [s.split('T')[0] for s in dates]
days_unique = np.unique(days)
days_unique = sorted(days_unique)
days_unique = days_unique[63:]

# print out csv title
print("Beginning calculation of values for csv file: " + csv_name)

# List of header strings
row_contents = ['date_obs', 'date_jd', 'v_quiet', 'v_disc', 'v_phot', 'v_conv', 'f_bright', 'f_spot', 'f', 'Bobs',
                'vphot_bright', 'vphot_spot', 'f_small', 'f_large', 'f_network', 'f_plage', 'f_nonconv',
                'quiet_flux', 'ar_flux', 'conv_flux', 'pol_flux', 'pol_conv_flux', 'vconv_quiet', 'vconv_large',
                'vconv_small']

# Append a list as new line to an old csv file
utils.append_list_as_row(csv_name, row_contents)

# get hmi data products
time_range = datetime.timedelta(seconds=22)
physobs_list = [a.Physobs.los_velocity, a.Physobs.los_magnetic_field, a.Physobs.intensity]

for i, date in enumerate(days_unique[41:48]):
    neid_use = np.isin(days, date)
    dtime = date[0:4] + "-" + date[4:6] + '-' + date[6:8] + 'T12:00:00'
    date_str, date_obj, date_jd = utils.get_dates(dtime)

    rv_in_file = []
    dvrms = []
    files = np.array(spec_fits_files)[neid_use]
    for f in files:
        main_header = fits.getheader(f)
        date_jd = fits.getheader(f)['OBSJD']

        if main_header['DRIFTFUN'] == 'simultaneous' and main_header['WAVECAL'] == 'LFCplusThAr':
            # check that uncertainty is not terrible
            if fits.getheader(f, 'CCFS')['DVRMS'] <= .0004:
                # get information from fits header
                rv_in_file.append(fits.getheader(f, 'CCFS')['CCFRVMOD'])
                dvrms.append(fits.getheader(f, 'CCFS')['DVRMS'])

    if len(rv_in_file) >= 3:
        # get the average (binned) rv value
        rv_bin = np.average(rv_in_file)
        rv_error = np.average(dvrms)

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
                # coordinate transformation for maps
                x, y, pd, r, d, mu = ctfuncs.coordinates(vmap)
                wij, nij, rij = ctfuncs.vel_coords(x, y, pd, r, vmap)

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
                # vrot = sfuncs.solar_rot_vel(wij, nij, rij, deltaw, deltan, deltar, dij, vmap, a_parameters)
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

                # create threshold array
                thresh_arr = sfuncs.thresh_map(fac_inds, spot_inds)

                # full threshold maps
                map_full_thresh = sfuncs.corrected_map(thresh_arr, mmap, map_type='Threshold',
                                                       frame=frames.HeliographicCarrington)

                ### velocity contribution due to convective motion of quiet-Sun
                v_quiet = sfuncs.v_quiet(map_vel_cor, imap, quiet)

                ### velocity contribution due to rotational Doppler imbalance of active regions (faculae/sunspots)
                # calculate photospheric velocity
                v_phot, vphot_bright, vphot_spot = sfuncs.v_phot(quiet, active, Lij, vrot, imap, mu, fac_inds, spot_inds)

                ### velocity contribution due to suppression of convective blueshift by active regions
                # calculate disc-averaged velocity
                v_disc = sfuncs.v_disc(map_vel_cor, imap)

                # calculate convective velocity
                v_conv = v_disc - v_quiet

                ### filling factor
                # calculate filling factor
                f_bright, f_spot, f = sfuncs.filling_factor(mu, mmap, active, fac_inds, spot_inds)

                ### unsigned magnetic flux
                # unsigned observed flux
                unsigned_obs_flux = sfuncs.unsigned_flux(map_mag_obs, imap)

                ### calculate the area filling factor
                pixA_hem = ctfuncs.pix_area_hem(wij, nij, rij, vmap)
                area = sfuncs.area_calc(active, pixA_hem)
                f_small, f_large, f_network, f_plage, f_nonconv = sfuncs.area_filling_factor(active, area, mu, mmap,
                                                                                             fac_inds)

                ### get the unsigned flux
                quiet_flux, ar_flux, conv_flux, pol_flux, pol_conv_flux = sfuncs.area_unsigned_flux(map_mag_obs, imap,
                                                                                                    area,
                                                                                                    active)

                ### get area weighted convective velocities
                vconv_quiet, vconv_large, vconv_small = sfuncs.area_vconv(map_vel_cor, imap, active, area)

                # make array of what we want to save
                save_vals = [v_quiet, v_disc, v_phot, v_conv, f_bright, f_spot, f, unsigned_obs_flux, vphot_bright,
                             vphot_spot, f_small, f_large, f_network, f_plage, f_nonconv, quiet_flux, ar_flux,
                             conv_flux, pol_flux, pol_conv_flux, vconv_quiet, vconv_large, vconv_small]

                # round stuff
                rounded = np.around(save_vals, 3)
                round_vals = [date_str, date_jd, rv_bin, rv_error]
                for val in rounded:
                    round_vals.append(val)

                # append these values to the csv file
                utils.append_list_as_row(csv_name, round_vals)
    else:
        # append these values to the csv file
        utils.append_list_as_row(bad_dates_file, 'Not enough RVs')
    # print that the date is completed
    print('\nCalculations and save to file complete for ' + date_str + ' index: ' + str(i))


# print elapsed time
end_time = time.time()

print((end_time - start_time) / 60)
