"""
Tamar Ervin
Date: July 1, 2021

determination of latitude based area threshold from
Milbourne 2019

"""

import numpy as np
import pandas as pd

from sklearn.linear_model import LinearRegression

import matplotlib.pyplot as plt

# csv file with rv components
csv_file = '/Users/tervin/csv_files/milbourne/area_fraction_cutoff_as_function_of_latitude.csv'

df = pd.read_csv(csv_file)

area = df.loguHem.values
sin_colat = df.sin_colat.values

# get best fit line
x = np.pi/2 - np.arcsin(sin_colat)
y = area
arr = np.polyfit(x, y, 2)
x_fit = np.linspace(x[0], x[-1], 100)
y_fit = arr[0]*x_fit**2 + arr[1]*x_fit + arr[2]

plt.scatter(x, y, label='Data', color='orchid')
plt.plot(x_fit, y_fit, label='Best Fit Line', color='lightblue')
# plt.plot(np.abs(lat), area_cuts, label='From data')
plt.xlabel('$\Theta$')
plt.ylabel("Area $\mu Hem$")
plt.title('Area cutoff as function of $\Theta$')
plt.legend()
plt.show()
plt.savefig('/Users/tervin/images/area_cutoff_sinlat_function.png')

#
# # get it to pixels
# r_sun = mmap.fits_header['rsun_obs']
# pix_dim = mmap.fits_header['cdelt1']
#
# sun_pix_area = np.pi*(r_sun)**2 / (pix_dim ** 2)
# hem_area = sun_pix_area/2 * 10e-6
#
# y_pix = y_fit * hem_area
#
# # mag thresh function -- get the map from save values
# from skimage.measure import label, regionprops
#
# Br_cutoff = 24
# mu_cutoff = 0.3
# mu = mmu
#
# # calculate magnetic cutoff value
# br_cutoff = Br_cutoff / mu
#
# # get active regions
# brthresh = np.full(shape=map_mag_cor.data.shape, fill_value=np.nan)
# use_indices = np.logical_and(mu > mu_cutoff, np.abs(map_mag_cor.data) > br_cutoff)
#
# # find isolated pixels
# # convert to integer array for skimage use
# y = use_indices.astype(int)
# # get area
# y_labeled = label(y, connectivity=2, background=0)
# y_area = [props.area for props in regionprops(y_labeled)]
#
# # get latitude
# mean_lat = np.rint([props.centroid for props in regionprops(y_labeled)]).astype(int)
# theta = np.arccos(mu)
# lat = [mtheta[l[0], l[1]] for l in mean_lat]
#
# # get areas corresponding to latitudes
# area_cuts = [(arr[0]*l**2 + arr[1]*l + arr[2])*hem_area for l in np.abs(lat)]
#
# # area constraint
# good_area = np.where(np.array(y_area) > area_cuts)
# good_area = good_area[0] + 1
# active_indices = np.isin(y_labeled, good_area)
#
# # add the relevant magnetic data to the thresholded array
# brthresh[active_indices] = map_mag_cor.data[active_indices]
#
# # create weights array
# quiet_indices = np.logical_and(mu > mu_cutoff, ~active_indices)
# mag_weights = np.full(shape=map_mag_cor.data.shape, fill_value=np.nan)
# mag_weights[quiet_indices] = 1
# mag_weights[active_indices] = 0
#
# plt.imshow(mag_weights)
# plt.colorbar()
# plt.show()