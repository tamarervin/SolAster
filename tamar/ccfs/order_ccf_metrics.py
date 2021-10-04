#!/usr/bin/env python
"""
Tamar Ervin
Date: August 26, 2021

Calculation of ccf metrics for individual orders.

Will be specifically looking at integrated area over
time as this seems to be the strongest proxy for magnetic
activity!!

"""

import sys

sys.path.append('/Users/tervin/NEID_Solar_analysis')

import os
import glob
import datetime
import pandas as pd
import numpy as np
from tqdm import tqdm
from scipy import constants
from sunpy.net import attrs as a

from astropy.io import fits

from scipy.integrate import simps

import NEIDcode

import tamar.tools.ccf_funcs as ccffuncs
from tamar.tools.settings import CsvDir, Config
import tamar.tools.utilities as utils

# file names
file_names = ['ccfs', 'above_weights', 'below_weights', '8weights' '7weights', '6weights', '78weights', '67weights', '56weights']
file_name = 'ccfs'

# get list of fits files
all_files = glob.glob(os.path.join(CsvDir.NEID_HOUR, '*.fits'))
dates = [d.split('/')[-1] for d in all_files]
dates = [d.split('_')[-1] for d in dates]
dates = [d.split('T')[0] for d in dates]

# read in SDO calculations csv
sdo_csv = os.path.join(CsvDir.NEID_CALC, 'rvs_from_fits.csv')
df = pd.read_csv(sdo_csv)
sdo_dates = df.date_obs.values

# mask to use
mask_name = file_name + '.txt'
if file_name == 'ccfs':
    mask_name = 'G2_espresso.txt'

# csv file to save ccf information
csv_name = os.path.join(CsvDir.CCFS, 'orders', file_name + '.csv')

# parameters for SDO/HMI image generation
time_range = datetime.timedelta(seconds=22)
physobs_list = [a.Physobs.los_velocity, a.Physobs.los_magnetic_field, a.Physobs.intensity]

# setup constants for CCF calculation
config = Config.config

# lightspeed
LIGHTSPEED = constants.c / 1000  # Speed of light in km/s
minzb = -30 / LIGHTSPEED
maxzb = +30 / LIGHTSPEED

# barycentric FWHM correction parameters -- KPNO
obsname = 'KPNO'
lat = 31.958092
longi = -111.600562
alt = 2091.0

# Set up velocity loop
velocity_loop = np.arange(config['velocity_min'], config['velocity_max'], config['velocity_step']) + config['qrv']
velsize = len(velocity_loop)
vel_span = 5

# generate mask object based on config parameters
mask = NEIDcode.Mask(mask_name, config['mask_half_width'], config['mask_environment'])
fsr_mask = CsvDir.FSR_MASK
fsr_mask = fits.getdata(fsr_mask)
fsr_mask = 1 - fsr_mask

# weights for building full CCF
ccf_weights = np.load(CsvDir.WEIGHTS)

# add row names to csv
row_contents = ['date_obs', 'date_jd', 'Bobs']
nums = [str(n) for n in np.arange(0, 122)]
row_contents += nums
utils.append_list_as_row(csv_name, row_contents)

# make lists for saving
for date in sdo_dates:
    if date == '2021-01-07T12:00:00':
        pass
    else:
        # get neid files
        date_use = date[0:4] + date[5:7] + date[8:10]
        date_use = np.isin(dates, date_use)
        files = np.array(all_files)[date_use]
        nfiles = len(files)
        print('Calculating Weighted CCF for', date, 'using', nfiles, 'files.')

        # total number of extracted orders
        orders = 122  # hard coded for NEID

        # make areas list
        areas = np.zeros(orders)

        # setup CCF calculation
        ccfs_fortran = np.zeros([nfiles, orders, velsize])
        ccfs_pipeline = np.zeros([nfiles, orders, velsize])
        rv_order = np.zeros([nfiles, orders])
        ccfs = np.zeros([orders, velsize])
        file_ccfs = np.zeros([nfiles, velsize])

        # calculate CCFs and RVs
        rv_in_file = []
        dvrms = []
        for n, f in enumerate(files):
            fd = fits.open(f)
            flux = fd['SCIFLUX'].data
            wave = fd['SCIWAVE'].data
            head = fd[0].header
            ccfhead = fd['CCFS'].header
            ccfwt = np.zeros(122)
            main_header = fits.getheader(f)
            date_jd = fits.getheader(f)['OBSJD']

            if main_header['DRIFTFUN'] == 'simultaneous' and main_header['WAVECAL'] == 'LFCplusThAr':
                # get RV and error from SDO calculations
                date_obs = date[0:4] + "-" + date[5:7].zfill(2) + '-' + date[8:10].zfill(2) + 'T12:00:00'
                sdo_inds = np.isin(sdo_dates, date_obs)
                rv = df.rv_model.values[sdo_inds] / 1e3
                B = df.Bobs.values[sdo_inds]
                vconv = df.v_conv.values[sdo_inds]
                jd_date = df.date_jd.values[sdo_inds]

                # check to even see if we should continue
                if len(rv) == 0:
                    pass
                else:
                    # check that uncertainty is not terrible
                    if fits.getheader(f, 'CCFS')['DVRMS'] <= .0004:
                        # switch between order index and echelle order, then loop
                        rv_in_file.append(fits.getheader(f, 'CCFS')['CCFRVMOD'])
                        dvrms.append(fits.getheader(f, 'CCFS')['DVRMS'])
                        for trueorder in tqdm(range(config['ordmin'], config['ordmax'] + 1, 1)):
                            zb = head['SSBZ%03i' % trueorder]
                            berv = head['SSBRV%03i' % trueorder]
                            raworder = config['bluest_order'] - trueorder
                            ccfwt[raworder] = ccfhead['CCFWT%03i' % trueorder]

                            # You have to remove NaNs ahead of passing spectrum to Fortran code
                            nanfree = NEIDcode.remove_nans(flux[raworder, :], method='linear')
                            spectrum = nanfree[0]

                            if fsr_mask is None:
                                # Number of blaze edge pixels to clip on either side of order
                                pix_start = config['clip_edge_pixels']
                                pix_end = np.shape(flux)[1] - config['clip_edge_pixels']
                            else:
                                if np.nansum(np.logical_not(fsr_mask[raworder, :])) == 0:
                                    print('fsr mask order summed to zero...sad')
                                    continue
                                else:
                                    pix_start = np.min(np.argwhere(np.logical_not(fsr_mask[raworder, :])))
                                    pix_end = np.max(np.argwhere(np.logical_not(fsr_mask[raworder, :])))

                            dummy_line_start = mask.start * ((1.0 + (velocity_loop[0] / LIGHTSPEED)) / (1.0 + maxzb))
                            dummy_line_end = mask.end * ((1.0 + (velocity_loop[-1] / LIGHTSPEED)) / (1.0 + minzb))
                            try:
                                line_index = np.where((dummy_line_start > np.min(wave[raworder, pix_start:pix_end])) &
                                                      (dummy_line_end < np.max(wave[raworder, pix_start:pix_end])))[0]
                            except TypeError:
                                print('Line index is None...devastating')
                                line_index = []

                            sn = np.ones(len(flux[raworder, pix_start:pix_end]))
                            if not len(line_index) == 0:
                                for k in range(len(velocity_loop)):
                                    ccfs_fortran[n, raworder, k] = NEIDcode.CCF_3d.ccf(mask.start[line_index],
                                                                                       mask.end[line_index],
                                                                                       wave[raworder, pix_start:pix_end],
                                                                                       spectrum[pix_start:pix_end],
                                                                                       mask.weight[line_index],
                                                                                       sn, velocity_loop[k], berv, 0.)
                            else:
                                ccfs_fortran[n, raworder, :] = np.zeros(len(velocity_loop))

        # average over all the files
        ccfs[:, :] = np.average(ccfs_fortran, axis=0)

        # get gaussian fit
        for j in range(orders):
            ccf = ccfs[j, :]
            ccf /= np.nansum(ccf)
            gaussian_fit, g_x, g_y, final_rv = ccffuncs.fit_gaussian_to_ccf(velocity_loop, ccf, config['qrv'],
                                                                            config['velocity_half_range_to_fit'])

            # shift ccf through interpolation
            ccf_shift = ccffuncs.shift_ccf(velocity_loop, ccf, final_rv)

            # get integration bounds
            low_bound, high_bound = np.argwhere(velocity_loop == - vel_span)[0][0], np.argwhere(velocity_loop == vel_span)[0][0]
            ccf_int = ccf_shift[low_bound:high_bound+1]
            int_area = simps(1 - ccf_int)
            areas[j] = int_area

        # append these values to the csv file
        vals = [date_obs, date_jd, B]
        for area in areas:
            vals.append(area)
        utils.append_list_as_row(csv_name, vals)

