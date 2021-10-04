#!/usr/bin/env python
"""
Tamar Ervin
Date: July 12, 2021

calculation of solar rvs and comparison to order
by order RVs
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

# name of csv file to get dates
dates_csv = os.path.join(CsvDir.NEID_CALC, 'rv_calcs.csv')

# get dates
data_df = pd.read_csv(dates_csv)
dates = data_df.date_obs.values

# get dates list
# path to NEID solar folder
neid_days = os.listdir(CsvDir.NEID_SOLAR)

# name of csv file to store calculations
csv_name = os.path.join(CsvDir.NEID_CALC, 'order_rvs.csv')

# name of csv to store bad dates
bad_dates_csv = os.path.join(CsvDir.NEID_BAD_DATES, 'long_bad_dates.csv')

# # # get pieces for date
# year = [d[0:4] for d in neid_days]
# month = [d[4:6] for d in neid_days]
# day = [d[6:8] for d in neid_days]
# # hour = [d[9:11] for d in days_list]
# # minute = [d[11:13] for d in days_list]
# # second = [d[-2:] for d in days_list]
# dates = [d[0:8] for d in neid_days]
# dates = np.unique(dates)

# get order numbers
order_numbers = np.arange(57, 170)
order_numbers = [str(n) for n in order_numbers]
order_numbers = [n.zfill(3) for n in order_numbers]

# print out csv title
print("Beginning calculation of values for csv file: " + csv_name)

# # List of header strings
obj_list = ['date_obs', 'date_jd', 'rv_sun', 'rv_error']
obj_list.reverse()
row_contents = order_numbers
for ele in obj_list:
    row_contents = [ele] + row_contents

# # Append a list as new line to an old csv file
utils.append_list_as_row(csv_name, row_contents)

# # get hmi data products
# time_range = datetime.timedelta(seconds=22)
# physobs_list = [a.Physobs.los_velocity, a.Physobs.los_magnetic_field, a.Physobs.intensity]

for i, date in enumerate(dates):

    file_date = date[0:4] + date[5:7] + date[8:10]
    file = os.path.join(CsvDir.NEID_SOLAR, file_date, 'level2', file_date)
    spec_fits_files = [i for i in glob.glob(os.path.join(file, '*.fits'))]
    # dtime = date[0:4] + "-" + date[4:6] + '-' + date[6:8] + 'T12:00:00.00'
    date_str, date_obj, date_jd = utils.get_dates(date)

    rv_in_file = []
    dv_in_file = []
    orders_in_file = []
    for f in spec_fits_files:

        # get date info
        date_jd = fits.getheader(f)['OBSJD']
        main_header = fits.getheader(f)

        # get RVs
        if main_header['DRIFTFUN'] == 'simultaneous' and main_header['WAVECAL'] == 'LFCplusThAr':
            # check that uncertainty is not terrible
            if fits.getheader(f, 'CCFS')['DVRMS'] <= .0004:
                # get information from fits header
                rv_in_file.append(fits.getheader(f, 'CCFS')['CCFRVMOD'])
                dv_in_file.append(fits.getheader(f, 'CCFS')['DVRMS'])
                order_rvs = [fits.getheader(f, 'CCFS')['CCFRV' + n] for n in order_numbers]
                orders_in_file.append(order_rvs)

    orders_use = list(zip(*orders_in_file))
    # get averages of order rvs, regular rvs, and uncertainty
    rv_bin = np.average(rv_in_file)
    dvrms = np.average(dv_in_file)
    order_rvs = [np.average(o) for o in orders_use]

    # append these values to the csv file
    save_list = [date_str, date_jd, rv_bin, dvrms]
    save_list.reverse()
    save_vals = order_rvs
    for ele in save_list:
        save_vals = [ele] + save_vals
    utils.append_list_as_row(csv_name, save_vals)

    # print completion statement
    print('\nCalculations and save to file complete for ' + date_str + ' index: ' + str(i))


# print elapsed time
end_time = time.time()

print((end_time - start_time) / 60)

