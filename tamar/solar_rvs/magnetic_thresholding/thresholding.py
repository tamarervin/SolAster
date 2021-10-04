"""
Tamar Ervin
Date: June 30, 2021

messing with thresholding for magnetic and intensity data

"""

import os

import sys

sys.path.append('//')

import time
import datetime
import numpy as np

import astropy.units as u

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.coordinates import frames

import tamar.tools.solar_funcs as sfuncs
import tamar.tools.coord_funcs as ctfuncs
import tamar.tools.lbc_funcs as lbfuncs
import tamar.tools.utilities as utils

# get hmi data products
cadence = a.Sample(24 * u.hour)  # querying cadence
start_date = datetime.datetime(2011, 1, 1, 12, 0, 0)
end_date = datetime.datetime(2011, 1, 2, 12, 0, 0)
time_range = datetime.timedelta(seconds=22)
cadence = a.Sample(24*u.hour)
physobs_list = [a.Physobs.los_velocity, a.Physobs.los_magnetic_field, a.Physobs.intensity]


# convert the date to a string -- required for use in csv file
date_str_start, date_obj_start = utils.get_dates(start_date)
date_str_end, date_obj_end = utils.get_dates(end_date)

# pull image within specified time range
result = Fido.search(a.Time(date_str_start, date_str_end),
                     a.Instrument.hmi, physobs_list[0] | physobs_list[1] | physobs_list[2], cadence)

# add file to list
file_download = Fido.fetch(result)
file_download1 = file_download[:2]
file_download1.append(file_download[-1])
# convert to map sequence
map_seq = sunpy.map.Map(sorted(file_download1))

# check for missing data types
missing_map = False
# split into data types
for i, map_obj in enumerate(map_seq):
    if map_obj.meta['content'] == 'DOPPLERGRAM':
        vmap = map_obj
    elif map_obj.meta['content'] == 'MAGNETOGRAM':
        mmap = map_obj
    elif map_obj.meta['content'] == 'CONTINUUM INTENSITY':
        imap = map_obj
    else:
        missing_map = True

# # transformation for magnetic maps
# mtheta, mphi, mmu = ctfuncs.coord_trans(mmap, R0=mmap.meta['rsun_obs'] * 760e3 / mmap.meta['rsun_ref'],
#                                         outside_map_val=np.nan)

mtheta, mphi, mmu = ctfuncs.coord_trans(mmap, R0=1.01,
                                        outside_map_val=np.nan)
good_mu = np.logical_and(mmu > 0.1, mmu != np.nan)
mmap.data[~good_mu] = np.nan
imap.data[~good_mu] = np.nan

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
Bobs, Br = sfuncs.mag_field(mmu, mmap, B_noise)

# corrected observed magnetic data map
map_mag_obs = sfuncs.corrected_map(Bobs, mmap, map_type='Corrected-Magnetogram',
                                   frame=frames.HeliographicCarrington)

# radial magnetic data map
map_mag_cor = sfuncs.corrected_map(Br, mmap, map_type='Corrected-Magnetogram',
                                   frame=frames.HeliographicCarrington)

# calculate magnetic threshold
# magnetic threshold value (G) from Yeo et. al. 2013
Br_cutoff = 24
# mu cutoff value
mu_cutoff = 0.3
# area cutoff value (pixels)
area_cutoff = 45
Brthresh, mag_weights = sfuncs.mag_thresh(mmu, map_mag_cor, Br_cutoff=Br_cutoff,
                                                mu_cutoff=mu_cutoff, area=area_cutoff)

# calculate intensity threshold
faculae, sunspots, fac_ind, spot_ind, int_cutoff = sfuncs.int_thresh(map_int_cor, map_mag_cor, mag_weights)

# create threshold array
thresh_arr = sfuncs.thresh_map(fac_ind, spot_ind)

# full threshold maps
map_full_thresh = sfuncs.corrected_map(thresh_arr, mmap, map_type='Threshold', frame=frames.HeliographicCarrington)

import matplotlib.pyplot as plt

plt.imshow(map_full_thresh.data, cmap=plt.get_cmap('bwr'))
plt.title("Comparison of faculae (blue) and sunspots (red) to magnetic data")
plt.savefig('/Users/tervin/images/area_mag_colors')

# area cutoff value (pixels)
area_cutoff = 1
Brthresh, mag_weights = sfuncs.mag_thresh(mmu, map_mag_cor, Br_cutoff=Br_cutoff,
                                                mu_cutoff=mu_cutoff, area=area_cutoff)

# calculate intensity threshold
faculae, sunspots, fac_ind, spot_ind, int_cutoff = sfuncs.int_thresh(map_int_cor, map_mag_cor, mag_weights)

# create threshold array
thresh_arr = sfuncs.thresh_map(fac_ind, spot_ind)

# full threshold maps
map_full_thresh = sfuncs.corrected_map(thresh_arr, mmap, map_type='Threshold', frame=frames.HeliographicCarrington)

plt.imshow(map_full_thresh.data, cmap=plt.get_cmap('bwr'))
plt.title("Comparison of faculae (blue) and sunspots (red) to magnetic data")
plt.savefig('/Users/tervin/images/no_area_mag_colors')


### cluster using kmeans
# from sklearn.cluster import KMeans
#
# X = np.abs(map_mag_cor.data)
# inds = np.logical_and(~np.isnan(X), mmu > 0)
# X = X[inds]
# arr = X.reshape(-1, 1)
# kmeans = KMeans(n_clusters=2, random_state=0, init='k-means++').fit(arr)
# labels = kmeans.labels_
#
# low_ind = np.where(labels == 1)
# high_ind = np.where(labels == 0)
# low_vals = X[low_ind]
# high_vals = X[high_ind]
#
# print(low_vals.max(), low_vals.min())
# print(high_vals.max(), high_vals.min())
#
# mag_data = np.full(shape=map_mag_cor.data.shape, fill_value=np.nan)
# mag_data[inds] = np.where(np.abs(map_mag_cor.data)[inds] > low_vals.min(), 1, 0)
#
# import matplotlib.pyplot as plt
# plt.imshow(mag_data)
# plt.colorbar()
# plt.show()
