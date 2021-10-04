"""
Tamar Ervin
Date: July 6, 2021

Mapping magnetic threshold maps and other relevant observables
over AIA images.
"""


import os

import sys

sys.path.append('/Users/tervin/NEID_Solar_analysis')

import time
import datetime
import pandas as pd
import numpy as np

import astropy.units as u

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.coordinates import frames

import tamar.tools.solar_funcs as sfuncs
import tamar.tools.lbc_funcs as lbfuncs
import tamar.tools.coord_funcs as ctfuncs
import tamar.tools.utilities as utils

start_time = time.time()

# get hmi data products
time_range = datetime.timedelta(seconds=22)
physobs_list = [a.Physobs.los_velocity, a.Physobs.los_magnetic_field, a.Physobs.intensity]

date = datetime.datetime(2011, 11, 10, 1, 29, 0)


# convert the date to a string -- required for use in csv file
date_str, date_obj = utils.get_dates(date)

# pull image within specified time range
result = Fido.search(a.Time(str(date_obj - time_range), str(date_obj + time_range)),
                     a.Instrument.hmi, physobs_list[0] | physobs_list[1] | physobs_list[2])
print(result)

# add file to list
file_download = Fido.fetch(result)

# remove unusable file types
good_files = []
for file in file_download:
    name, extension = os.path.splitext(file)
    if extension == '.fits':
        good_files.append(file)

if len(good_files) != 3:
    # add the data
    # append these values to the csv file
    save_vals = [date_str, 'not three good files']

    pass
else:
    # convert to map sequence
    map_seq = sunpy.map.Map(sorted(good_files))

    # check for missing data types
    missing_map = False
    # split into data types
    for j, map_obj in enumerate(map_seq):
        if map_obj.meta['content'] == 'DOPPLERGRAM':
            vmap = map_obj
        elif map_obj.meta['content'] == 'MAGNETOGRAM':
            mmap = map_obj
        elif map_obj.meta['content'] == 'CONTINUUM INTENSITY':
            imap = map_obj
        else:
            missing_map = True

    if missing_map:
        print("Missing a data product for " + date_str)

        pass

    else:
        # transformation for velocity maps
        vtheta, vphi, vmu = ctfuncs.coord_trans(vmap, R0=vmap.meta['rsun_obs'] * 760e3 / vmap.meta['rsun_ref'], outside_map_val=np.nan)
        vmu, vmap = ctfuncs.fix_mu(vmu, vmap)

        # transformation for intensity maps
        itheta, iphi, imu = ctfuncs.coord_trans(imap, R0=imap.meta['rsun_obs'] * 760e3 / imap.meta['rsun_ref'], outside_map_val=np.nan)
        imu, imap = ctfuncs.fix_mu(imu, imap)

        # transformation for magnetic maps
        mtheta, mphi, mmu = ctfuncs.coord_trans(mmap, R0=mmap.meta['rsun_obs'] * 760e3 / mmap.meta['rsun_ref'], outside_map_val=np.nan)
        mmu, mmap = ctfuncs.fix_mu(mmu, mmap)

        # calculate relative positions
        deltaw, deltan, deltar, dij = sfuncs.rel_positions(vtheta, vphi, vmap)

        # calculate spacecraft velocity
        vsc = sfuncs.spacecraft_vel(deltaw, deltan, deltar, dij, vmap)

        # optimized solar rotation parameters
        a1 = 14.713
        a2 = -2.396
        a3 = -1.787
        a_parameters = [a1, a2, a3]

        # calculation of solar rotation velocity
        vrot = sfuncs.solar_rot_vel(a_parameters, vtheta, vphi, deltaw, deltan, deltar, dij, vmap)

        # calculate corrected velocity
        corrected_vel = vmap.data - vrot - vsc

        # corrected velocity maps
        map_vel_cor = sfuncs.corrected_map(corrected_vel, vmap, map_type='Corrected-Dopplergram',
                                           frame=frames.HeliographicCarrington)

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

        ### calculate magnetic threshold
        # magnetic threshold value (G) from Yeo et. al. 2013
        Br_cutoff = 24
        # mu cutoff value
        mu_cutoff = 0.3

        # calculate magnetic threshold
        Brthresh, mag_weights = sfuncs.mag_thresh(mmu, map_mag_cor, mmap, Br_cutoff=Br_cutoff, mu_cutoff=mu_cutoff)

        # create array of active weights (flipped from mag weights)
        act_weights = sfuncs.active_weights(mag_weights)

        # thresholded magnetic maps
        map_mag_thresh = sfuncs.corrected_map(Brthresh, mmap, map_type='Magnetic-Threshold',
                                              frame=frames.HeliographicCarrington)

        # calculate intensity threshold
        faculae, sunspots, fac_ind, spot_ind, int_cutoff = sfuncs.int_thresh(map_int_cor, map_mag_cor, mag_weights)

        # create threshold array
        thresh_arr = sfuncs.thresh_map(fac_ind, spot_ind)

        # full threshold maps
        map_full_thresh = sfuncs.corrected_map(thresh_arr, mmap, map_type='Threshold', frame=frames.HeliographicCarrington)

        ### velocity contribution due to convective motion of quiet-Sun
        v_quiet = sfuncs.v_quiet(map_vel_cor, imap, mag_weights)

        ### velocity contribution due to rotational Doppler imbalance of active regions (faculae/sunspots)
        # calculate photospheric velocity
        v_phot = sfuncs.v_phot(mag_weights, act_weights, Lij, vsc, imap, vmap)

        ### velocity contribution due to suppression of convective blueshift by active regions
        # calculate disc-averaged velocity
        v_disc = sfuncs.v_disc(map_vel_cor, imap)

        # calculate convective velocity
        v_conv = v_disc - v_quiet - v_phot

        ### filling factor
        # calculate filling factor
        filling_factor = sfuncs.filling_factor(map_mag_obs, act_weights)

        # calculate sunspot and faculae filling factors
        f_bright, f_spot = sfuncs.int_filling_factor(fac_ind, spot_ind)

        ### unsigned magnetic flux
        # unsigned observed flux
        unsigned_obs_flux = sfuncs.unsigned_flux(map_mag_obs, imap)

        # make array of what we want to save
        save_vals = [v_quiet,  v_disc, v_phot, v_conv, f_bright, f_spot, filling_factor, unsigned_obs_flux]

        # round stuff
        rounded = np.around(save_vals, 3)
        round_vals = [str(mmap.meta['date-obs'])]
        for val in rounded:
            round_vals.append(val)

        # print that the date is completed
        print('\nCalculations complete for ' + date_str)

#### ----- get aia image
# pull image within specified time range
time_range = datetime.timedelta(seconds=6)
result = Fido.search(a.Time(str(date_obj - time_range), str(date_obj + time_range)),
                     a.Instrument.aia, a.Wavelength(1600*u.Angstrom))
file_download = Fido.fetch(result)
aia = sunpy.map.Map(sorted(file_download))

#### ----- plot the stuff
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# plot intensity map
good_int = np.where(imu < 0.3, np.nan, map_int_cor.data)
int_map = sfuncs.corrected_map(good_int, imap, map_type='Corrected-Intensitygram',
                                           frame=frames.HeliographicCarrington)
int_map.plot_settings['cmap'] = plt.get_cmap('hinodesotintensity')
int_map.plot_settings['norm'] = colors.Normalize()
int_map.plot()
plt.colorbar()
plt.show()

# plot velocity map
good_vel = np.where(vmu < 0.3, np.nan, vmap.data)
vel_map = sfuncs.corrected_map(good_vel, vmap, map_type='Corrected-Dopplergram',
                                   frame=frames.HeliographicCarrington)
vel_map.plot_settings['cmap'] = plt.get_cmap('Greys')
vel_map.plot()
plt.colorbar()
plt.show()

# plot magnetic map
good_mag = np.where(mmu < 0.3, np.nan, map_mag_cor.data)
mag_map = sfuncs.corrected_map(good_mag, mmap, map_type='Corrected-Magnetogram',
                                   frame=frames.HeliographicCarrington)
mag_map.plot_settings['cmap'] = plt.get_cmap('hmimag')
mag_map.plot_settings['norm'] = colors.Normalize()
mag_map.plot()
plt.colorbar()
plt.show()

# plot threshold map
good_thresh = np.where(mmu < 0.35, np.nan, map_full_thresh.data)
thresh_map = sfuncs.corrected_map(good_thresh, mmap, map_type='Threshold-Map',
                                   frame=frames.HeliographicCarrington)
mag_map.plot_settings['cmap'] = plt.get_cmap('hmimag')
mag_map.plot_settings['norm'] = colors.Normalize()
mag_map.plot()
thresh_map.plot_settings['cmap'] = plt.get_cmap('bwr_r')
thresh_map.plot()
plt.colorbar()
plt.show()

#### making subplots
fig, axs = plt.subplots(2, 2)
fig.suptitle("November 10, 2011 at 1:30 UTC")

# fix data
good_int = np.where(imu < 0.3, np.nan, map_int_cor.data)
good_mag = np.where(mmu < 0.3, np.nan, map_mag_cor.data)
good_vel = np.where(mmu < 0.3, np.nan, vmap.data)
good_thresh = np.where(mmu < 0.3, np.nan, map_full_thresh.data)

# make maps
vel_map = sfuncs.corrected_map(good_vel, vmap, map_type='Corrected-Dopplergram',
                                   frame=frames.HeliographicCarrington)
int_map = sfuncs.corrected_map(good_int, imap, map_type='Corrected-Intensitygram',
                                           frame=frames.HeliographicCarrington)
mag_map = sfuncs.corrected_map(good_mag, mmap, map_type='Corrected-Magnetogram',
                                   frame=frames.HeliographicCarrington)
thresh_map = sfuncs.corrected_map(good_thresh, mmap, map_type='Threshold-Map',
                                   frame=frames.HeliographicCarrington)

# rotate hmi data
map_int_cor_rot = int_map.rotate(order=3)
map_mag_cor_rot = mag_map.rotate(order=3)
vmap_rot = vel_map.rotate(order=3)
map_full_thresh_rot = thresh_map.rotate(order=3)


# plot data
axs[0, 0].imshow(map_int_cor_rot.data, cmap=plt.get_cmap('hinodesotintensity'), norm=colors.Normalize())
axs[0, 1].imshow(vmap_rot.data, cmap=plt.get_cmap('Greys'),  norm=colors.Normalize())
axs[1, 0].imshow(map_mag_cor_rot.data, cmap=plt.get_cmap('hmimag'), norm=colors.Normalize())
axs[1, 1].imshow(map_full_thresh_rot.data, cmap=plt.get_cmap('bwr_r'))

# titles
axs[0, 0].set_title("Flattened Intensity")
axs[0, 1].set_title("Original Velocity Map")
axs[1, 0].set_title("Unsigned Flux")
axs[1, 1].set_title("Threshold Map")

# remove axes
for j in range(0, 2):
    for k in range(0, 2):
        axs[j, k].set_xticks([])
        axs[j, k].set_yticks([])

# # colorbars -- TODO: need scalar mappable colormap?
# fig.colorbar(ax=axs[0, 0])
# fig.colorbar(ax=axs[0, 1])
# fig.colorbar(ax=axs[1, 0])

plt.show()
# plt.savefig('/Users/tervin/images/12_30_2020/hmi_plot')


### overlay magnetic threshold map
plt.title('Magnetic Threshold')
plt.imshow(map_mag_cor_rot.data, cmap=plt.get_cmap('hmimag'), norm=colors.Normalize())
plt.colorbar()
plt.imshow(map_full_thresh_rot.data, cmap=plt.get_cmap('bwr_r'))
plt.show()
# plt.savefig('/Users/tervin/images/12_30_2020/mag_thresh_overlay')

### overlay threshold map on the actual sun!

# # add file to list
#
# plt.title('Solar Magnetic Thresholding')
aia.plot()
plt.show()
# map_thresh = map_full_thresh.rotate(order=3)
# use_thresh = sfuncs.corrected_map(map_thresh.data, mmap, map_type='Threshold-Map',
#                                    frame=frames.Helioprojective)
# plt.imshow(aia.data, cmap=plt.get_cmap('sdoaia1600'), norm=colors.LogNorm())
# plt.imshow(use_thresh.data, cmap='bwr_r')
# # plt.imshow(aia[0].data, cmap=plt.get_cmap('sdoaia211'), norm=colors.LogNorm())
# plt.show()
