"""
Tamar Ervin
Date: July 8, 2021

updates to coordinate transformations to work with
solar rotational velocity calculation

"""

import numpy as np
import datetime
from sunpy.net import attrs as a
import sunpy.map
from sunpy.net import Fido

# read in map
time_range = datetime.timedelta(seconds=22)
date_obj = datetime.datetime(2021, 1, 1, 12, 00, 00)
# pull image within specified time range
result = Fido.search(a.Time(str(date_obj - time_range), str(date_obj + time_range)),
                     a.Instrument.hmi, a.Physobs.los_velocity)

# add file to list
file_download = Fido.fetch(result)

# get map
smap = sunpy.map.Map(file_download)

# parameters
outside_map_val = np.nan
obsv_lon = smap.fits_header['crln_obs']
obsv_lat = smap.fits_header['crlt_obs']
R0=smap.meta['rsun_obs'] * 760e3 / smap.meta['rsun_ref']


# get x and y coordinates
header = smap.fits_header

# pixel array
axis1 = np.linspace(1, header['NAXIS1'], header['NAXIS1']) - header['CRPIX1']
axis2 = np.linspace(1, header['NAXIS2'], header['NAXIS2']) - header['CRPIX2']

# pixel offset from center of image, x,y
x, y = np.meshgrid(axis1, axis2)

# distance of observer to sun in solar radii
d = smap.fits_header['DSUN_OBS'] / smap.fits_header['RSUN_REF']

# focal length in pixels
f = 180. * 3600. / np.pi / smap.fits_header['CDELT1']

# distance (in pixels) to pixel
pd = np.sqrt(x ** 2 + y ** 2)

# distance (in solar r) to pixel
pr = f * f * pd * pd + pd ** 4 - d * d * pd ** 4 + 0.J
r = (d * f * pd - np.sqrt(pr)) / (f * f + pd * pd)

# separate complex parts
r = r.real

# get mu array
pr = 1 - r ** 2 + 0.J
cos_alpha = (np.sqrt(pr)).real
sin_alpha = r.real
cos_theta = ((d - np.sqrt(pr)) / np.sqrt(r ** 2 + (d - np.sqrt(pr)) ** 2)).real
sin_theta = (np.sqrt(1 - cos_theta ** 2)).real
mu = cos_alpha * cos_theta - sin_alpha * sin_theta

# get the west, north, and radial velocity arrays
head = smap.fits_header
crota2 = head['CROTA2']  # deg

# transform each pixel to get into Heliographic CR Coordinates
dw = y * np.sin(np.deg2rad(crota2)) + x * np.cos(np.deg2rad(crota2))
dn = y * np.cos(np.deg2rad(crota2)) - x * np.sin(np.deg2rad(crota2))

# get cartesian coordinates for velocity calculations
pr = 1 - r ** 2 + 0.J
wij = r * dw / pd
nij = r * dn / pd
rij = (np.sqrt(pr)).real