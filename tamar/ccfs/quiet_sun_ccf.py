"""
Tamar Ervin
Date: July 19, 2021

Creation of quiet-Sun template
CCF for each order
- quiet-Sun criteria:
    - good weather
    - no sunspots (low f total and f spot = 0)
    - low flux (total, active)
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits

from tamar.tools.settings import CsvDir, ImgDir

# ignore division error
np.seterr(divide='ignore', invalid='ignore')

# read in magnetic calculations csv
csv_name = 'rvs_from_fits.csv'
csv_file = os.path.join(CsvDir.NEID_CALC, csv_name)
df = pd.read_csv(csv_file)

# read in dates and magnetic information
dates = df.date_obs.values
f = df.f.values
f_bright = df.f_bright.values
f_spot = df.f_spot.values
Bobs = df.Bobs.values
flux_ar = df.ar_flux.values

# first cut: no sunspots
no_spots = np.where(f_spot == 0)

print('Date with lowest magnetic flux:', dates[np.argmin(Bobs)])
print('Date with lowest filling factor:', dates[np.argmin(f)])
print('Dates with no sunspots:', dates[no_spots])

# second cut: "low" total filling factor -- one sigma below median
print('The max filling factor (%) for days with no sunspots:', np.nanmax(f[no_spots]))
print('The min filling factor (%) for days with no sunspots:', np.nanmin(f[no_spots]))
low_f = np.where(f <= np.nanmedian(f[no_spots]) - np.std(f[no_spots]))
print("Dates with low filling factor and no sunspots:", dates[low_f])

# third cut: magnetic flux
print('The max unsigned flux (G) for days with no sunspots:', np.nanmax(Bobs[no_spots]))
print('The min unsigned flux (G) for days with no sunspots:', np.nanmin(Bobs[no_spots]))
print('Median Unsigned Flux (G):', np.nanmedian(Bobs))
low_flux = np.where(Bobs <= np.nanmedian(Bobs[low_f]) - np.std(Bobs[low_f]))
print("Dates with low flux, filling factor, and no sunspots:", dates[low_flux])

# fourth cut: lowest active flux
low_ar = np.argmin(flux_ar[low_flux])
quiet_day = dates[low_flux[0][low_ar]]
print('The quiet Sun day used will be:', quiet_day)

# get list of directories to build structure
neid_days = os.listdir(CsvDir.NEID_SOLAR)

# get date of quiet-Sun template
date = quiet_day[0:4] + quiet_day[5:7] + quiet_day[8:10]

# get paths to files
spec_path = os.path.join(CsvDir.NEID_SOLAR, date, 'level2', date)

# search for all fits files within 'spec_path'
spec_fits_files = [i for i in glob.glob(os.path.join(spec_path, '*.fits'))]
if len(spec_fits_files) < 1:
    print('no fits files found!')

# get CCFs for all files
file_ccfs = []
for f in spec_fits_files:

    # get all CCFs
    ccfs = fits.getdata(f, 'CCFS')
    nord = ccfs.shape[0]

    # normalize order CCFs
    norm_orders = []
    for ind, order in enumerate(range(nord)):
        ccf = ccfs[ind]
        ccf_scaled = (ccf - np.amin(ccf))
        ccf_scaled /= max(ccf_scaled)
        norm_orders.append(ccf_scaled)

    # create CCFs list for files
    file_ccfs.append(norm_orders)

# average CCFs for each order
averaged_orders = []
for ind, order in enumerate(range(nord)):
    order_ccfs = [f[ind] for f in file_ccfs]
    averaged = np.nanmean(order_ccfs, axis=0)
    averaged_orders.append(averaged)

#### ----- PLOTTING ----- ####
# set plot parameters
Col = plt.get_cmap('Spectral')
plt.style.use('seaborn-darkgrid')
plt.rcParams['figure.figsize'] = [5, 7.5]

# get the velocity information from fits file
vel_start = fits.getheader(spec_fits_files[0], 'CCFS')['CCFSTART']
vel_step = fits.getheader(spec_fits_files[0], 'CCFS')['CCFSTEP']
nvels = fits.getheader(spec_fits_files[0], 'CCFS')['NAXIS1']
vel_arr = np.arange(vel_start, vel_start + nvels * vel_step, vel_step)

# get wavelength information from fits file
fib = 'SCI'
fits_extension_wavelength = fib + 'WAVE'
wvl = fits.getdata(spec_fits_files[0], fits_extension_wavelength) / 10.

# plot each order cross-correlation function
for ind, order in enumerate(range(nord)):
    wvl_norm = 1. - (np.mean(wvl[ind,:]) - 420.) / (720. - 420.)
    ccf = averaged_orders[ind]
    ccf_scaled = (ccf - np.amin(ccf))
    ccf_scaled /= max(ccf_scaled) * 0.3
    ccf_scaled += ind
    plt.plot(vel_arr, ccf_scaled, c='black', linewidth=3)
    plt.plot(vel_arr, ccf_scaled, c=Col(wvl_norm), lw=2.)
    plt.annotate(str(173 - ind), (max(vel_arr) + 2., ccf_scaled[-1] - 0.5), fontsize=7)

plt.title('Normalized Quiet Sun Order CCFs:\n' + quiet_day)
plt.xlabel('Velocity [km s$^{-1}$]')
plt.show()
plt.savefig(os.path.join(ImgDir.CCF_IMG, 'quiet_sun_ccf_' + date + '.png'))