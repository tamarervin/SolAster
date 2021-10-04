"""
Tamar Ervin
Date: June 21, 2021

calculation of unsigned flux
"""

import numpy as np

import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord


import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.coordinates import frames

import tamar.tools.solar_funcs as sfuncs

# get magnetograms
cadence = a.Sample(24*u.hour)  # querying cadence
start_date = '2021-01-01T12:00:00'  # start date of query
end_date = '2021-01-05T12:00:00'
physobs_list = [a.Physobs.los_velocity, a.Physobs.los_magnetic_field, a.Physobs.intensity]

result = Fido.search(a.Time(start_date, end_date),
                     a.Instrument.hmi, physobs_list[0] | physobs_list[1] | physobs_list[2], cadence)


# download results
file_download = Fido.fetch(result)
# convert to map sequence
map_seq = sunpy.map.Map(sorted(file_download))
hmi_vel, hmi_mag, hmi_int \
    = [], [], []
for i, map_obj in enumerate(map_seq):
    if map_obj.meta['content'] == 'DOPPLERGRAM':
        hmi_vel.append(map_obj)
    elif map_obj.meta['content'] == 'MAGNETOGRAM':
        hmi_mag.append(map_obj)
    elif map_obj.meta['content'] == 'CONTINUUM INTENSITY':
        hmi_int.append(map_obj)

# choose one map for starters
mag_map = hmi_mag[0]

# calculate mu value (from solar_funcs.py)
mag_map_rot = mag_map.rotate(order=3)
x, y = sfuncs.get_scales_from_map(mag_map_rot)
R0 = mag_map_rot.meta['rsun_obs'] * 760e3 / mag_map_rot.meta['rsun_ref']
theta, phi, mu = sfuncs.get_coordinates(x, y, mag_map_rot, R0=R0, outside_map_val=np.nan)

# magnetic data
Bobs = mag_map_rot.data
# get valid indices
use_indices = np.logical_and(mu > 0, mu != np.nan)
# calculate full magnetic field strength
Br = np.full(shape=mag_map_rot.data.shape, fill_value=np.nan)
Br[use_indices] = Bobs[use_indices] / mu[use_indices]
# TODO: fix calculation of mu value, def not working right

# make corrected map
coord = SkyCoord(0*u.arcsec, 0*u.arcsec, obstime=mag_map_rot.date, observer=mag_map_rot.observer_coordinate, frame=frames.HeliographicCarrington)
header = sunpy.map.make_fitswcs_header(Br, coord)
header['content'] = mag_map_rot.meta['content']
header['telescop'] = mag_map_rot.meta['telescop']
header['wavelnth'] = mag_map_rot.meta['wavelnth']
corrected_map = sunpy.map.Map(Br, header)
corrected_map.peek()
mag_map_rot.peek()

# magnetic thresholding
Br_cutoff = 24/mu
mu_cutoff = 0.3
Brthresh = np.full(shape=mag_map_rot.data.shape, fill_value=np.nan)
use_indices = np.logical_and(mu > mu_cutoff, Br > Br_cutoff)
Brthresh[use_indices] = Br[use_indices]
plt.imshow(Brthresh)
plt.show()

# get weights array
mag_weight = np.ones(shape=mag_map_rot.data.shape)  # TODO: do we need to account for outside of limb pixels
mag_weight[use_indices] = 0
plt.imshow(mag_weight)
plt.show()

# intensity thresholding
int_map = hmi_int[0].rotate(order=3)
### theoretically we would correct this for limb brightening
# for now we will just use the regular intensity
Iflat = int_map.data
# calculate quiet sun intensity
int_quiet = np.nansum(Iflat*mag_weight) / np.nansum(mag_weight)
# intensity threshold
int_cutoff = 0.89 * int_quiet

# threhsold by intensity
faculae = np.full(shape=int_map.data.shape, fill_value=np.nan)
sunspots = np.full(shape=int_map.data.shape, fill_value=np.nan)
fac_ind = np.logical_and(Br > Br_cutoff, Iflat > int_cutoff)
spot_ind = np.logical_and(Br > Br_cutoff, Iflat < int_cutoff)

faculae[fac_ind] = Br[fac_ind]
sunspots[spot_ind] = Br[spot_ind]

plt.imshow(faculae)
plt.show()

plt.imshow(sunspots)
plt.show()

