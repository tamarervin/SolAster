"""
Tamar Ervin
Date: June 21, 2021

Correction the magnetogram data for foreshortening and
calculating the unsigned magnetic field.

"""

import numpy as np

import matplotlib.pyplot as plt

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.coordinates import frames

from skimage.measure import label, regionprops

import sdo_hmi_rvs.tools.calculation_funcs as sfuncs
import sdo_hmi_rvs.tools.coord_funcs as ctfuncs

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


# coordinate transformation for maps
x, y, pd, r, d, mu = ctfuncs.coordinates(vmap)
wij, nij, rij = ctfuncs.vel_coords(x, y, pd, r, vmap)

# remove bad mu values
vmap, mmap, imap = ctfuncs.fix_mu(mu, [vmap, mmap, imap])
# remove bad mu values
vmap, mmap, imap = ctfuncs.fix_mu(mu, [vmap, mmap, imap])

# calculate relative positions
deltaw, deltan, deltar, dij = sfuncs.rel_positions(wij, nij, rij, vmap)

# magnetic noise level
B_noise = 8
mu_cutoff = 0.3

# get valid indices
use_indices = np.logical_and(mu > mu_cutoff, mu != np.nan)
mag_indices = np.logical_and(use_indices, np.abs(mmap.data) < B_noise)

# calculate full magnetic field strength
Bobs = mmap.data
Br = np.full(shape=mmap.data.shape, fill_value=np.nan)
Br[use_indices] = Bobs[use_indices] / mu[use_indices]
Bobs[mag_indices] = 0
Br[mag_indices] = 0
