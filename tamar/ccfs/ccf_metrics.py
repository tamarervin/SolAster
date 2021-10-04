#!/usr/bin/env python
"""
Tamar Ervin
Date: August 2, 2021

calculation of RVs and metrics from CCFs
using Arpita's function for Gaussian fit

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

# file name
file_name = 'espresso+neid_mask_97_to_108'

# for file_name in file_names:
# mask to use
mask_name = file_name + '.mas'
if file_name == 'ccfs':
    mask_name = 'G2_espresso.txt'

# csv file to save ccf information
pickle_csv = os.path.join(CsvDir.CCFS, file_name + '.csv')

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
vel_span = 5

# generate mask object based on config parameters
mask = NEIDcode.Mask(mask_name, config['mask_half_width'], config['mask_environment'])
fsr_mask = CsvDir.FSR_MASK
fsr_mask = fits.getdata(fsr_mask)
fsr_mask = 1 - fsr_mask

# get list of fits files
all_files = glob.glob(os.path.join(CsvDir.NEID_HOUR, '*.fits'))
dates = [d.split('/')[-1] for d in all_files]
dates = [d.split('_')[-1] for d in dates]
dates = [d.split('T')[0] for d in dates]

# read in SDO calculations csv
sdo_csv = os.path.join(CsvDir.NEID_CALC, 'rvs_from_fits.csv')
df = pd.read_csv(sdo_csv)
sdo_dates = df.date_obs.values

# make lists for saving
csv_dates, jd_dates, rv_sun, rv_model, error, ccfs, shift_ccfs, gaussian, rv_gauss, skew, area, \
Bobs, v_conv = [], [], [], [], [], [], [], [], [], [], [], [], []
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

        # setup CCF calculation
        ccfs_fortran = np.zeros([nfiles, orders, velsize])
        ccfs_pipeline = np.zeros([nfiles, orders, velsize])
        rv_order = np.zeros([nfiles, orders])
        ccf_weighted = np.zeros([nfiles, velsize])
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
                print('Date is good:', date)
                # get RV and error from SDO calculations
                date_obs = date[0:4] + "-" + date[4:6].zfill(2) + '-' + date[6:8].zfill(2) + 'T12:00:00'
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

                # check again that we have sdo data
                if len(rv) == 0:
                    pass
                else:
                    # calculate full CCF
                    ccfsum = np.nansum(ccfs_fortran[n, :, :], axis=0)

                    # apply weighting scheme to each order to get ccf from each file
                    ccfmod = np.nansum((ccfs_fortran[n, :, :].T * ccfwt).T, axis=0)
                    rescale = np.nansum(ccfsum) / np.nansum(ccfmod)
                    ccfmod *= rescale
                    file_ccfs[n, :] = ccfmod

        # sum up all of the weighted CCFs to make the 'final' CCF
        ccf = np.nansum(file_ccfs, axis=0)

        # normalize
        ccf /= np.nanmax(ccf)

        # get gaussian fit
        gaussian_fit, g_x, g_y, final_rv = ccffuncs.fit_gaussian_to_ccf(velocity_loop, ccf, config['qrv'],
                                                                        config['velocity_half_range_to_fit'])

        # shift ccf through interpolation
        ccf_shift = ccffuncs.shift_ccf(velocity_loop, ccf, final_rv)

        # get integration bounds
        low_bound, high_bound = np.argwhere(velocity_loop == - vel_span)[0][0], np.argwhere(velocity_loop == vel_span)[0][0]
        ccf_int = ccf_shift[low_bound:high_bound+1]
        int_area = simps(1 - ccf_int)

        # get integrated area of ccf
        # int_area, int_error = ccffuncs.integrated_area(gaussian_fit)

        # calculate skew for each CCF -- skew = 3 * (mean - median)/std
        skw = 3*(np.mean(ccf) - np.median(ccf)) / np.std(ccf)

        # add these to lists
        csv_dates.append(date)
        jd_dates.append(jd_date)
        rv_model.append(rv)
        ccfs.append(ccf)
        gaussian.append(gaussian_fit)
        rv_gauss.append(final_rv)
        Bobs.append(B)
        v_conv.append(vconv)
        rv_sun.append(np.mean(rv_in_file))
        error.append(np.mean(dvrms))
        skew.append(skw)
        area.append(int_area)

        # print calculation complete statement
        print('Calculation complete for', date)

# create dataframe
d = {
    'dates': csv_dates,
    'jd_dates': jd_dates,
    'rv_model': rv_model,
    'rv_error': error,
    'rv_sun': rv_sun,
    'ccf': ccfs,
    'gaussian': gaussian,
    'rv_gauss': rv_gauss,
    'median_skew': skew,
    'integrated_area': area,
    'Bobs': Bobs,
    'v_conv': v_conv
}

# save dataframe to csv using pickle
pickle_df = pd.DataFrame(data=d)
pickle_df.to_pickle(pickle_csv)

# ## get neid ccf
# rvs = []
# for n, f in enumerate(files):
#     ccfs_neid = np.zeros([nfiles, orders, velsize])
#     vel_start = fits.getheader(f, 'CCFS')['CCFSTART']
#     vel_step = fits.getheader(f, 'CCFS')['CCFSTEP']
#     nvels = fits.getheader(f, 'CCFS')['NAXIS1']
#     vel_arr = np.arange(vel_start, vel_start + nvels * vel_step, vel_step)
#     ccfs_neid[n, :, :] = fits.getdata(f, 'CCFS')
#     rvs.append(fits.getheader(f, 'CCFS')['CCFRVMOD'])
#     for ord_n in range(122):
#         ccf_sum = np.nansum(ccfs_neid[n, ord_n, :])
#         ccf_sums.append(ccf_sum)
#         # apply weighting scheme to each order
#         ccfs_scaled[ord_n, :] = (ccfs_neid[n, ord_n, :] / ccf_sum) * ccf_weights[ord_n]
#
#     # sum up all of the weighted CCFs to make the 'final' CCF
#     ccf_weighted[n, :] = np.nansum(ccfs_scaled, axis=0)
#
# # daily binned CCF
# neid_ccf = np.average(ccf_weighted, axis=0)
# # normalize
# neid_ccf /= np.nanmax(neid_ccf)
# # neid RVs
# neid_rv = np.mean(rvs)
#
# # plot the stuff
# from astropy.modeling import models, fitting
# import matplotlib.pyplot as plt
#
# plt.plot(velocity_loop, ccf, label='Fortran', color='lightpink')
# p = models.Gaussian1D(amplitude=gaussian_fit.amplitude, mean=gaussian_fit.mean, stddev=gaussian_fit.stddev)
# plt.plot(velocity_loop, p(velocity_loop), color='plum', label='Fit')
# plt.legend()
# # plt.title('Fit Compared to CCF\nNEID RV: ' + str(np.round(neid_rv, 3)) +
# #           ' Fit RV: ' + str(np.round(gaussian_fit.mean.value, 3)) +
# #           ' SDO RV: ' + str(np.round(rv[0], 3)))
# plt.show()
#
# # percentage off for calculated RV
# diff = gaussian_fit.mean.value - neid_rv
# print('Error:', np.abs(diff) / neid_rv * 100)
#
#
# # looking at fit ccf
# # gaussian = A exp[-(x - mean)^2 / (2*std^2)]
#
#
# def gauss(x, mean, std, A):
#     g = A * np.exp(-(x - mean) ** 2 / (2 * std ** 2))
#     return g
#
#
# x = velocity_loop
# A = gaussian_fit.amplitude
# mean = gaussian_fit.mean
# std = gaussian_fit.stddev
# g = gauss(x, mean, std, A)
# plt.plot(velocity_loop, g)
# plt.show()
