#!/usr/bin/env python
"""
Tamar Ervin
Date: August 5, 2021

calculation of RVs and metrics from CCFs
using Arpita's function for Gaussian fit
for April 17, 2021

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

import NEIDcode

import tamar.tools.ccf_funcs as ccffuncs
from tamar.tools.settings import CsvDir, Config

# mask to use
mask_name = 'G2_espresso.txt'

# csv file to save ccf information
pickle_csv = os.path.join(CsvDir.CCFS, '04_17_ccfs.csv')

# get list of dates!
dates = os.listdir(CsvDir.NEID_SOLAR)
dates = sorted(dates)

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

# generate mask object based on config parameters
mask = NEIDcode.Mask(mask_name, config['mask_half_width'], config['mask_environment'])
fsr_mask = '/Users/tervin/NEID_Solar_analysis/NEIDcode/masks/neidMaster_FSR_Mask20210331_v001.fits'
fsr_mask = fits.getdata(fsr_mask)
fsr_mask = 1 - fsr_mask

# make lists for saving
rv_sun, rv_error, ccfs, gaussian, rv_gauss, jd_dates = [], [], [], [], [], []

# get neid files
files = [i for i in glob.glob(os.path.join('/Users/tervin/solar_cme_day_04_17', '*.fits'))]
nfiles = len(files)
print('Calculating Weighted CCF for solar CME day April 14, 2021 using', nfiles, 'files.')

# total number of extracted orders
orders = 122  # hard coded for NEID

# setup CCF calculation
ccfs_fortran = np.zeros([nfiles, orders, velsize])
ccfs_pipeline = np.zeros([nfiles, orders, velsize])
rv_order = np.zeros([nfiles, orders])
ccf_weighted = np.zeros([nfiles, velsize])

# calculate CCFs and RVs
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
        # check that uncertainty is not terrible
        if fits.getheader(f, 'CCFS')['DVRMS'] <= .0004:
            # get information from fits header
            rv_sun.append(fits.getheader(f, 'CCFS')['CCFRVMOD'])
            rv_error.append(fits.getheader(f, 'CCFS')['DVRMS'])
        # check to even see if we should continue
        if len(rv_sun) == 0:
            pass
        else:
            # check that uncertainty is not terrible
            if fits.getheader(f, 'CCFS')['DVRMS'] <= .0004:
                # switch between order index and echelle order, then loop
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

            # calculate full CCF
            ccfs_scaled = np.zeros_like(ccfs_fortran[n, :, :])
            ccfsum = np.nansum(ccfs_fortran[n, :, :], axis=0)

            # apply weighting scheme to each order
            ccfmod = np.nansum((ccfs_fortran[n, :, :].T * ccfwt).T, axis=0)
            rescale = np.nansum(ccfsum) / np.nansum(ccfmod)
            ccfmod *= rescale

            # get gaussian fit
            gaussian_fit, g_x, g_y, final_rv = ccffuncs.fit_gaussian_to_ccf(velocity_loop, ccfmod,
                                                                            config['qrv'],
                                                                            config['velocity_half_range_to_fit'])
            print('Gaussian fit for file number', n)
        # add these to lists
        ccfs.append(ccfmod)
        gaussian.append(gaussian_fit)
        rv_gauss.append(final_rv)
        jd_dates.append(date_jd)

# create dataframe
d = {
    'dates': jd_dates,
    'rv_sun': rv_sun,
    'rv_error': rv_error,
    'ccf': ccfs,
    'gaussian': gaussian,
    'rv_gauss': rv_gauss
}

# save dataframe to csv using pickle
pickle_df = pd.DataFrame(data=d)
pickle_df.to_pickle(pickle_csv)
