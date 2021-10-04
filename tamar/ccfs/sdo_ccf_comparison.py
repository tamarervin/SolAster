"""
Tamar Ervin
Date: July 23, 2021

- looking at CCFs in comparison to SDO/HMI images
on a specific date
- comparing residual CCFs to these observations
"""

import os
import glob
import datetime
import pandas as pd
from tqdm import tqdm
from scipy import constants

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import NEIDcode

import tamar.tools.utilities as utils
import tamar.tools.ccf_funcs as ccffuncs
from tamar.tools.settings import CsvDir, ImgDir, PlotDir, Config

from astropy.coordinates import EarthLocation
from astropy.time import Time
from barycorrpy.PhysicalConstants import *

# make list of dates of interest
dates = ['20210205', '20201226', '20210225', '20210327', '20210423']

# csv file to save ccf information
pickle_csv = os.path.join(CsvDir.CCFS, 'ccf_pickle.csv')

# parameters for SDO/HMI image generation
time_range = datetime.timedelta(seconds=22)
physobs_list = [a.Physobs.los_velocity, a.Physobs.los_magnetic_field, a.Physobs.intensity]

# setup constants for CCF calculation
config = Config.config

# barycentric FWHM correction parameters -- KPNO
obsname = 'KPNO'
lat = 31.958092
longi = -111.600562
alt = 2091.0

# lightspeed
LIGHTSPEED = constants.c / 1000  # Speed of light in km/s
minzb = -30 / LIGHTSPEED
maxzb = +30 / LIGHTSPEED

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

    # setup CCF calculation
    # total number of extracted orders
    orders = 122  # hard coded for NEID
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
                # ccfs_fortran[n, :, :] = (ccfs_fortran[n, :, :] ** 2. - delta) ** 0.5

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

            # get RV and error from SDO calculations
            date_obs = date[0:4] + "-" + date[4:6].zfill(2) + '-' + date[6:8].zfill(2) + 'T12:00:00'
            sdo_inds = np.isin(sdo_dates, date_obs)
            rv_model = df.rv_model_weather_coeff.values[sdo_inds] / 1e3
            rv_error = df.rv_error.values[sdo_inds] / 1e3

            # shift CCF
            ccf_weighted[n, :] = ccffuncs.shift_ccf(velocity_loop, ccf_weighted[n, :], rv_model)

            # add the collier cameron correction
            # fwhm /= (fwhm - delta)
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

df = pd.DataFrame(data=d)
# save dataframe to csv using pickle
df.to_pickle(pickle_csv)

# plot these things!!

# convert date to sunpy usable form
dtime = [dt[0:4] + "-" + dt[4:6] + '-' + dt[6:8] + 'T12:01:00' for dt in dates[1:]]
dts = [utils.get_dates(dt) for dt in dtime]
date_obj = [dt[1] for dt in dts]
date_str = [dt[0] for dt in dts]

# pull image within specified time range
result = [Fido.search(a.Time(str(dt - time_range), str(dt + time_range)),
                      a.Instrument.hmi, physobs_list[0] | physobs_list[1] | physobs_list[2]) for dt in date_obj]

# add file to list
file_download = [Fido.fetch(r) for r in result]

# make map sequence
map_seq = sunpy.map.Map(file_download)

# split into data types
vmap, mmap, imap = [], [], []
for j, map_obj in enumerate(map_seq):
    if map_obj.meta['content'] == 'DOPPLERGRAM':
        vmap.append(map_obj)
    elif map_obj.meta['content'] == 'MAGNETOGRAM':
        mmap.append(map_obj)
    elif map_obj.meta['content'] == 'CONTINUUM INTENSITY':
        imap.append(map_obj)
    else:
        missing_map = True

# plot style
plot_style = os.path.join(PlotDir.MPL, 'timeseries.mplstyle')
plt.style.use(plot_style)
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.cal'] = 'Helvetica'

# make plot structure
fig, axs = plt.subplots(4, 4, figsize=[24, 24], gridspec_kw={'hspace': 0.2, 'wspace': 0.4})
fig.suptitle('SDO Images and Corresponding CCFs', fontsize=18, y=0.9)

# get CCFs and RVs from pickle file
fpickle = pd.read_pickle(pickle_csv)
ccf_list = fpickle.ccf.values
rv = fpickle.rv_model.values

# get quiet-Sun CCF
quiet_ccf = ccf_list[0]

# get quiet-Sun RV
quiet_rv = rv[0]

# plot stuff!
for i in range(len(ccf_list) - 1):
    # get date
    use_date = date_str[i].split('T')

    # get SDO RV
    residual_rv = rv[i + 1] - quiet_rv

    # start by adding the intensity plots - row one
    axs[i, 0].imshow(imap[i].data, cmap=plt.get_cmap('hinodesotintensity'), norm=colors.Normalize())
    axs[i, 0].set_title('Intensitygram: ' + str(use_date[0]))
    axs[i, 0].set_xticks([])
    axs[i, 0].set_yticks([])
    axs[i, 0].set_box_aspect(1)
    axs[i, 0].set_facecolor('snow')
    axs[i, 0].grid()

    # then add the magnetogram plots - row two
    mmap[i].plot_settings['norm'] = plt.Normalize(-100, 100)
    axs[i, 1].imshow(mmap[i].data, cmap=plt.get_cmap('hmimag'))
    axs[i, 1].set_title('Magnetogram: ' + str(use_date[0]))
    axs[i, 1].set_xticks([])
    axs[i, 1].set_yticks([])
    axs[i, 1].set_box_aspect(1)
    axs[i, 1].set_facecolor('snow')
    axs[i, 1].grid()

    # add the CCFs - row three
    axs[i, 2].plot(velocity_loop, ccf_list[i + 1], c='black', linewidth=3)
    axs[i, 2].set_title('CCF: ' + str(use_date[0]))
    axs[i, 2].set_xlabel('Velocity [km s$^{-1}$]')
    axs[i, 2].set_ylabel('CCF')
    axs[i, 2].set_xlim(-25, 25)
    axs[i, 2].set_box_aspect(1)
    axs[i, 2].set_facecolor('snow')

    # calculate and add residual CCFs - row four
    residual = ccf_list[i + 1] - quiet_ccf
    axs[i, 3].plot(velocity_loop, residual, c='black', linewidth=3)
    axs[i, 3].set_title('Residual CCF: ' + str(use_date[0]))
    axs[i, 3].set_xlabel('Velocity [km s$^{-1}$]')
    axs[i, 3].set_ylabel('Residual CCF')
    axs[i, 3].set_xlim(-25, 25)
    axs[i, 3].set_ylim(min(residual) - 0.00005, max(residual) + 0.00005)
    axs[i, 3].set_box_aspect(1)
    axs[i, 3].set_facecolor('snow')

    # add residual RV information
    text = 'SDO Derived Residual RV: ' + str(np.round(residual_rv[0] * 1000, 3)) + ' m/s'
    props = dict(boxstyle='round', facecolor='lightpink', alpha=0.5)
    axs[i, 3].text(0.05, 0.95, text, transform=axs[i, 3].transAxes, fontsize=10,
                   verticalalignment='top', bbox=props)
    print('Residual Weighted Mean:', np.average((ccf_list[i + 1] - quiet_ccf), weights=velocity_loop))

# show and save plot
plt.show()
# save_fig = os.path.join(ImgDir.CCF_IMG, 'sdo_ccf_comparison.png')
# plt.savefig(save_fig)



