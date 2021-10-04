"""
Tamar Ervin
Date: July 22, 2021

- calculation of order by order RVs using
the Fortran code from Arpita/Sam
- not entirely sure of the best way to save
these calculations...
-
"""

import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import constants
from pathlib import Path

from tqdm import tqdm
from mpl_toolkits.mplot3d import Axes3D
import time
import scipy.interpolate as interp

from tamar.tools.settings import CsvDir

import NEIDcode

# constant parameters
config = {
    'velocity_min': -100,  # km/s minumum velocity at which to compute the CCF
    'velocity_max': 101,  # km/s  minumum velocity at which to compute the CCF
    'velocity_step': 0.25,  # km/s velocity step size at which to compute the CCF
    'qrv': -33.165,  # km/s systemic velocity of target star
    'mask_name': 'G2_espresso.txt',  # name of mask file containing line list
    'mask_half_width': 0.5,  # km/s width of binary mask (same for all lines)
    'mask_environment': 'air',  # air or vac wavelengths for mask lines
    'ordmin': 57,  # min echelle order number for which to compute CCF (reddest)
    'ordmax': 169,  # max echelle order number for which to compute CCF (bluest)
    'bluest_order': 173,  # bluest order currently extracted (don't change this)
    'resolution': 120000,  # spectral resolution of instrument, lambda / delta_lambda
    'sampling': 5,  # approximate pixel sampling of resolution element (pixels, FWHM)
    'velocity_half_range_to_fit': 20,
    # half width of velocity window for which to fit CCF and derive RV (subset of fill velocity range)
    'clip_edge_pixels': 500  # toss out pixels towards the edge of the array
}

# constants
LIGHTSPEED = constants.c / 1000  # Speed of light in km/s
minzb = -30 / LIGHTSPEED
maxzb = +30 / LIGHTSPEED

# Set up velocity loop
velocity_loop = np.arange(config['velocity_min'], config['velocity_max'], config['velocity_step']) + config['qrv']
velsize = len(velocity_loop)

# generate mask object based on config parameters
mask = NEIDcode.Mask(config['mask_name'], config['mask_half_width'], config['mask_environment'])

# fsr mask -- what is this...
fsr_mask = None

# get list of neid days
neid_days = os.listdir(CsvDir.NEID_SOLAR)

# calculate individual order CCFs for each day
ccf_days = []
for date in neid_days[:1]:

    # get files list
    file = os.path.join(CsvDir.NEID_SOLAR, date, 'level2', date)
    files = [i for i in glob.glob(os.path.join(file, '*.fits'))]
    nfiles = len(files)
    print('There were {} files found for processing for {}'.format(nfiles, date))

    # total number of extracted orders
    orders = 122  # hard coded for NEID
    ccfs_fortran = np.zeros([nfiles, orders, velsize])
    ccfs_pipeline = np.zeros([nfiles, orders, velsize])
    rv_order = np.zeros([nfiles, orders])

    start_time = time.time()

    for n, f in enumerate(files):
        a = fits.open(f)
        flux = a['SCIFLUX'].data
        wave = a['SCIWAVE'].data
        head = a[0].header
        ccfhead = a['CCFS'].header
        ccfwt = np.zeros(122)

        # switch between order index and echelle order, then loop
        for trueorder in tqdm(range(config['ordmin'], config['ordmax'] + 1, 1)):
            zb = head['SSBZ%03i' % trueorder]
            berv = head['SSBRV%03i' % trueorder]
            raworder = config['bluest_order'] - trueorder
            ccfwt[raworder] = ccfhead['CCFWT%03i' % trueorder]

            # remove NaNs ahead of passing spectrum to Fortran code
            nanfree = NEIDcode.remove_nans(flux[raworder, :], method='linear')
            spectrum = nanfree[0]

            if fsr_mask is None:
                # number of blaze edge pixels to clip on either side of order
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
                    ccfs_fortran[n, raworder, k] = NEIDcode.CCF_3d.ccf(mask.start[line_index], mask.end[line_index],
                                                                       wave[raworder, pix_start:pix_end],
                                                                       spectrum[pix_start:pix_end],
                                                                       mask.weight[line_index],
                                                                       sn, velocity_loop[k],
                                                                       berv, 0.)

            else:
                ccfs_fortran[n, raworder, :] = np.zeros(len(velocity_loop))

    ccf_days.append(ccfs_fortran)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print("Finished processing {} files in {} seconds.".format(nfiles, elapsed_time))
