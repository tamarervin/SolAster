"""
Tamar Ervin
Date: June 21, 2021

flatten continuum intensity images to correct for limb-darkening
use fifth order polynomial and constants from Allen 1973

IDL code here: https://hesperia.gsfc.nasa.gov/ssw/gen/idl/solar/darklimb_correct.pro

"""
import os
import sys

import astropy.units as u

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../'))
import SolAster.tools.coord_funcs as ctfuncs

# get magnetograms
cadence = a.Sample(24 * u.hour)  # querying cadence
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

imap = hmi_int[0].rotate(order=3)

# calculate mu value (from org_vel_corrections.py)
x, y = ctfuncs.get_scales_from_map(imap)
theta, phi, mu = ctfuncs.get_coordinates(x, y, imap, R0=1.01, outside_map_val=0)


# calculate limb-brightening correction
# limb brightening coefficients
def get_u(ll):
    """
    function to get limb darkening coefficient u (Allen 1973) based on wavelength value

    Parameters
    ----------
    ll: wavelength

    Returns
    -------
    vl: coefficient v based on wavelength
    """

    pll = np.array([1.0, ll, ll ** 2, ll ** 3, ll ** 4, ll ** 5])
    au = -8.9829751
    bu = 0.0069093916
    cu = -1.8144591e-6
    du = 2.2540875e-10
    eu = -1.3389747e-14
    fu = 3.0453572e-19
    a = np.array([au, bu, cu, du, eu, fu])
    ul = sum(a * pll)
    return ul


def get_v(ll):
    """
    function to get limb darkening coefficient v (Allen 1973) based on wavelength value

    Parameters
    ----------
    ll: wavelength

    Returns
    -------
    vl: coefficient v based on wavelength

    """

    pll = np.array([1.0, ll, ll ** 2, ll ** 3, ll ** 4, ll ** 5])
    av = 9.2891180
    bv = -0.0062212632
    cv = 1.5788029e-6
    dv = -1.9359644e-10
    ev = 1.1444469e-14
    fv = -2.599494e-19
    a = np.array([av, bv, cv, dv, ev, fv])
    vl = sum(a * pll)
    return vl


# get data
data = imap.data

# get header information
header = imap.fits_header

# wavelength
wavelength = header['WAVELNTH']

# coordinates
x_center = header['CRPIX1']
y_center = header['CRPIX2']
radius = header['RSUN_OBS'] / header['CDELT1']

# size
size = imap.data.shape[0]

# wavelength as float
ll = 1.0 * wavelength

# convert data to array and set nan values to zero
array = np.array(data)
NaNs = np.isnan(array)
array[NaNs] = 0.0

# get correction coefficients based on wavelength
ul = get_u(ll)
vl = get_v(ll)

# make arrays of x and y coordinates
x = np.arange(0, size, 1.0)
y = np.arange(0, size, 1.0)

# create grids for x, y arrays
x_mat, y_mat = np.meshgrid(x, y)

# make z array such that we get a circle
z = np.sqrt((x_mat - x_center) ** 2 + (y_mat - y_center) ** 2)

# normalize circle to radius of 1
grid = z / radius
out = np.where(grid > 1.0)
grid[out] = 0.0

# calculate polynomial
Lij = 1.0 - ul - vl + ul * np.cos(np.arcsin(grid)) + vl * np.cos(np.arcsin(grid)) ** 2

# flatten intensity array
Iflat = imap.data / Lij
