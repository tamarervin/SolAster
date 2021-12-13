"""
Tamar Ervin
Date: July 7, 2021

velocity corrections

"""

import os
import sys

import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.coordinates import frames

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../'))
import SolAster.tools.coord_funcs as ctfuncs

# get hmi data products
cadence = a.Sample(24 * u.hour)  # querying cadence
start_date = '2021-01-01T12:00:00'  # start date of query
end_date = '2021-01-03T12:00:00'
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
vmap = hmi_vel[0]
mmap = hmi_mag[0]
imap = hmi_int[0]

# transformation for velocity maps
# coordinate transformation for maps
px, py, p, r, d, mu = ctfuncs.coordinates(vmap)
rw_obs, rn_obs, rr_obs = ctfuncs.vel_coords(px, py, p, r, vmap)

# spacecraft velocity correction
# coordinates of each pixel
wij = rw_obs
nij = rn_obs
rij = rr_obs

# calculate relative positions of each pixel
rsc = vmap.meta['dsun_obs'] / vmap.meta['rsun_ref']
deltaw = wij
deltan = nij
deltar = rij - rsc
dij = np.sqrt(deltaw ** 2 + deltan ** 2 + deltar ** 2)

# velocity of spacecraft relative to sun
vscw = vmap.meta['obs_vw']
vscn = vmap.meta['obs_vn']
vscr = vmap.meta['obs_vr']

# pixel-wise magnitude of spacecraft velocity
vsc = - (deltaw * vscw + deltan * vscn + deltar * vscr) / dij

# solar rotation velocity correction

# apply to cartesian coordinates
x1 = wij
y1 = nij * np.cos(np.deg2rad(vmap.meta['crlt_obs'])) + rij * np.sin(np.deg2rad(vmap.meta['crlt_obs']))
z1 = -nij * np.sin(np.deg2rad(vmap.meta['crlt_obs'])) + rij * np.cos(np.deg2rad(vmap.meta['crlt_obs']))

hx = x1 * np.cos(np.deg2rad(vmap.meta['crln_obs'])) + z1 * np.sin(np.deg2rad(vmap.meta['crln_obs']))
hy = y1
hz = -x1*np.sin(np.deg2rad(vmap.meta['crln_obs'])) + z1*np.cos(np.deg2rad(vmap.meta['crln_obs']))

# differential velocity profile
a_parameters = [14.713, -2.396, -1.787]
w = (a_parameters[0] + a_parameters[1] * ((np.sin(hy)) ** 2) + a_parameters[2] * ((np.sin(hy)) ** 4)) * 1./86400.* np.pi/180.

# apply in rotation frame
vx_rot = w * hz * vmap.meta['rsun_ref']
vy_rot = 0.
vz_rot = -w * hx * vmap.meta['rsun_ref']

v1 = np.cos(np.deg2rad(vmap.meta['crln_obs']))*vx_rot - np.sin(np.deg2rad(vmap.meta['crln_obs']))*vz_rot
v2 = vy_rot
v3 = np.sin(np.deg2rad(vmap.meta['crln_obs']))*vx_rot + np.cos(np.deg2rad(vmap.meta['crln_obs']))*vz_rot

# project into correct direction
vrotw = v1
vrotn = v2*np.cos(np.deg2rad(vmap.meta['crlt_obs'])) - v3*np.sin(np.deg2rad(vmap.meta['crlt_obs']))
vrotr = v2*np.sin(np.deg2rad(vmap.meta['crlt_obs'])) + v3*np.cos(np.deg2rad(vmap.meta['crlt_obs']))

# get full rotational velocity
vrot = - (deltaw * vrotw + deltan * vrotn + deltar * vrotr) / dij

# subtract from the Doppler images
corrected_vel = vmap.data - vrot - vsc

# make corrected map
coord = SkyCoord(0*u.arcsec, 0*u.arcsec, obstime=vmap.date, observer=vmap.observer_coordinate, frame=frames.HeliographicCarrington)
header = sunpy.map.make_fitswcs_header(corrected_vel, coord)
header['content'] = 'CORRECTED_DOPPLERGRAM'
header['telescop'] = vmap.meta['telescop']
header['wavelnth'] = vmap.meta['wavelnth']
corrected_map = sunpy.map.Map(corrected_vel, header)

# plot maps to compare
vmap.peek()
corrected_map.peek()