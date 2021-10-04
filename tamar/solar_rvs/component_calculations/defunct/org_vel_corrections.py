"""
Tamar Ervin
Date: June 21, 2021


outline for step 2 in calculation of solar RVs
based on method from Haywood et. al. 2016
- calculation of spacecraft motion velocity
"""

import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.coordinates import frames

import tamar.tools.solar_funcs as sfuncs

# get hmi data products
cadence = a.Sample(24*u.hour)  # querying cadence
start_date = '2021-01-01T12:00:00'  # start date of query
end_date = '2021-01-05T12:00:00'
physobs_list = [a.Physobs.los_velocity, a.Physobs.los_magnetic_field, a.Physobs.intensity]

result = Fido.search(a.Time(start_date, end_date),
                     a.Instrument.hmi, physobs_list[0] | physobs_list[1] | physobs_list[2], cadence)

# print results
# print(result)

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
vel_map = hmi_vel[0]

# rsun = np.float32(vel_map.rsun_meters)
rsun = vel_map.meta['rsun_obs'] * 760e3 / vel_map.meta['rsun_ref']
x, y = sfuncs.get_scales_from_map(vel_map)
theta, phi, mu = sfuncs.get_coordinates(x, y, vel_map, R0=rsun, outside_map_val=np.nan)

# coordinates of each pixel
wij = theta
nij = phi
rij = rsun

### spacecraft velocity
# calculate relative positions of each pixel
rsc = vel_map.meta['dsun_obs']
deltaw = wij
deltan = nij
deltar = rsun - rsc  # TODO: this gives a negative radial position...bad?? -- I think not
dij = np.sqrt(deltaw**2 + deltan**2 + deltar**2)

# velocity of spacecraft relative to sun
vscw = vel_map.meta['obs_vw']
vscn = vel_map.meta['obs_vn']
vscr = vel_map.meta['obs_vr']

# pixel-wise magnitude of spacecraft velocity
vsc = - (deltaw * vscw + deltan * vscn + deltar * vscr)/dij

### solar rotation
a1 = 14.713  # TODO: are there newer optimized parameters?? -- don't think there are...
a2 = -2.396
a3 = -1.787
# apply parameters to determine vrot for given image pixel
w = a1 + a2 * (np.sin(theta))**2 + a3 * (np.sin(theta))**4

# TODO: not sure if the application of this rotation is correct -- I think it is
vrotw = vscw * w
vrotn = vscn * w
vrotr = vscr * w
vrot = - (deltaw * vrotw + deltan * vrotn + deltar * vrotr)/dij

# subtract from the Doppler images
corrected_vel = vel_map.data - vrot - vsc

# make corrected map
coord = SkyCoord(0*u.arcsec, 0*u.arcsec, obstime=vel_map.date, observer=vel_map.observer_coordinate, frame=frames.HeliographicCarrington)
header = sunpy.map.make_fitswcs_header(corrected_vel, coord)
header['content'] = 'CORRECTED_DOPPLERGRAM'
header['telescop'] = vel_map.meta['telescop']
header['wavelnth'] = vel_map.meta['wavelnth']
corrected_map = sunpy.map.Map(corrected_vel, header)

# plot maps to compare
vel_map.peek()
corrected_map.peek()




