#!/usr/bin/env python
"""
Tamar Ervin
Date: July 19, 2021

Calculation of residual CCFs
for individual orders
- averaged over 30 minutes around
solar noon to remove effects of
p-mode oscillations
- differential CCF is the quiet-Sun
template for individual orders subtracted
from the averaged order CCFs

"""

import sys

sys.path.append('/Users/tervin/NEID_Solar_analysis')

import os
import glob
import datetime
import pandas as pd
from tqdm import tqdm
from scipy import constants
from sunpy.net import attrs as a

from astropy.io import fits

import NEIDcode

import tamar.tools.ccf_funcs as ccffuncs
from tamar.tools.settings import CsvDir, Config

from astropy.coordinates import EarthLocation
from astropy.time import Time
from barycorrpy.PhysicalConstants import *

# get list of dates!
dates = os.listdir(CsvDir.NEID_SOLAR)
dates = sorted(dates)

# csv file to save ccf information
pickle_csv = os.path.join(CsvDir.CCFS, 'ccf_pickle.csv')

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
mask = NEIDcode.Mask(config['mask_name'], config['mask_half_width'], config['mask_environment'])
fsr_mask = None

# weights for building full CCF
ccf_weights = np.load(CsvDir.WEIGHTS)

# read in SDO calculations csv
sdo_csv = os.path.join(CsvDir.NEID_CALC, 'rvs_from_txt.csv')
df = pd.read_csv(sdo_csv)
sdo_dates = df.date_obs.values

# make lists for saving
rv_sun, error, ccfs, shift_ccfs = [], [], [], []
for date in dates:

    # get neid file
    file = os.path.join(CsvDir.NEID_SOLAR, date, 'level2', date)
    files = [i for i in glob.glob(os.path.join(file, '*.fits'))]
    nfiles = len(files)
    print('Calculating Weighted CCF for', date, 'using', nfiles, 'files.')

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
            # get RV and error from SDO calculations
            # check to even see if we should continue
            date_obs = date[0:4] + "-" + date[4:6].zfill(2) + '-' + date[6:8].zfill(2) + 'T12:00:00'
            sdo_inds = np.isin(sdo_dates, date_obs)
            rv_model = df.rv_model_weather_coeff.values[sdo_inds] / 1e3
            rv_error = df.rv_error.values[sdo_inds] / 1e3
            if len(rv_model) == 0:
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
                            if np.sum(np.logical_not(fsr_mask[raworder, :])) == 0:
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

                    # account for barycentric correction
                    if obsname:
                        loc = EarthLocation.of_site(obsname)
                        lat = loc.lat.value
                        longi = loc.lon.value
                        alt = loc.height.value
                    else:
                        loc = EarthLocation.from_geodetic(longi, lat, height=alt)

                    JDOBJ = Time(date_jd, format='jd', scale='utc')
                    delta = np.array(ccffuncs.CalculateFWHMDifference_SolarRotation_Ecliptic(loc, JDOBJ))
                    ccfs_fortran[n, :, :] = (ccfs_fortran[n, :, :] ** 2. - delta) ** 0.5

            # check again that we have sdo data
            if len(rv_model) == 0:
                pass
            else:
                # calculate full CCF
                ccfs_scaled = np.zeros_like(ccfs_fortran[n, :, :])
                ccf_sums = []

                # loop through each order and normalize the CCFs to their integrated areas
                for ord_n in range(122):
                    ccf_sum = np.nansum(ccfs_fortran[n, ord_n, :])
                    ccf_sums.append(ccf_sum)
                    # apply weighting scheme to each order
                    ccfs_scaled[ord_n, :] = (ccfs_fortran[n, ord_n, :] / ccf_sum) * ccf_weights[ord_n]

                # sum up all of the weighted CCFs to make the 'final' CCF
                ccf_weighted[n, :] = np.nansum(ccfs_scaled, axis=0)

                # shift CCF
                ccf_weighted[n, :] = ccffuncs.shift_ccf(velocity_loop, ccf_weighted[n, :], rv_model)

    # daily binned CCF
    ccf = np.average(ccf_weighted, axis=0)
    # normalize
    ccf /= np.nanmax(ccf)

    # add these to lists
    rv_sun.append(rv_model)
    error.append(rv_error)
    ccfs.append(ccf)

    # print calculation complete statement
    print('Calculation complete for', date)

# create dataframe
d = {
    'dates': dates,
    'rv_model': rv_sun,
    'rv_error': error,
    'ccf': ccfs
}

# save dataframe to csv using pickle
pickle_df = pd.DataFrame(data=d)
pickle_df.to_pickle(pickle_csv)

# read in data
# get CCFs and RVs from pickle file
pickle_csv = os.path.join(CsvDir.CCFS, 'ccfs.csv')
fpickle = pd.read_pickle(pickle_csv)
dates = fpickle.dates.values
ccf_list = fpickle.ccf.values
rv = fpickle.rv_model.values
rv_error = fpickle.rv_error.values

# get quiet-Sun ccf
quiet = np.where(dates == '20210205')
quiet_ccf = ccf_list[quiet][0]
rv_quiet = rv[quiet]
ccf_list = np.delete(ccf_list, quiet)
rv = np.delete(rv, quiet)
rv_error = np.delete(rv_error, quiet)
dates = np.delete(dates, quiet)

# find bad CCFs -- need to figure out why this is...
bad_ccfs = np.array([np.isnan(ccf[0]) for ccf in ccf_list])
ccf_list = ccf_list[~bad_ccfs]
rv = rv[~bad_ccfs]
rv_error = rv_error[~bad_ccfs]
dates = dates[~bad_ccfs]
# setup constants for CCF calculation
config = Config.config

# calculate residuals and save to pickle file
residual_ccfs = [ccf - quiet_ccf for ccf in ccf_list]
rv_residual = [r - rv_quiet for r in rv]
residual_pickle = os.path.join(CsvDir.CCFS, 'residual_ccfs.csv')
d = {
    'dates': dates,
    'rv_residual': rv_residual,
    'rv_error': rv_error,
    'residual_ccf': residual_ccfs
}

# save dataframe to csv using pickle
pickle_df = pd.DataFrame(data=d)
pickle_df.to_pickle(residual_pickle)

