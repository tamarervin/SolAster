"""
Tamar Ervin
Date: June 23, 2021

Limb-darkening correction functions:
flatten continuum intensity images to correct for limb-darkening
use fifth order polynomial and constants from Allen 1973

"""

import numpy as np


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


def limb_polynomial(imap):
    """
    function to calculate limb darkening correction polynomial based on IDl function:
    https://hesperia.gsfc.nasa.gov/ssw/gen/idl/solar/darklimb_correct.pro

    Parameters
    ----------
    imap: UNCORRECTED Sunpy map object (Intensitygram)

    Returns
    -------
    Lij: limb-darkening polynomial array

    """

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

    return Lij


### defunct

# def get_dict(dict_name):
#     """
#     function to return the correct dictionary of coefficient values
#
#     Parameters
#     ----------
#     dict_name: name of the coefficient dictionary you want
#
#     Returns
#     -------
#     dict: coefficient dictionary of choosing
#
#     """
#
#     if dict_name == 'u2_coeff':
#         return u2_coeff
#     elif dict_name == 'v2_coeff':
#         return v2_coeff
#     else:
#         return limb_bright_coeff
#
#
# limb_bright_coeff = {
#     0.20: [0.12, 0.33],
#     0.22: [-1.3, 1.6],
#     0.245: [-0.1, 0.85],
#     0.265: [-0.1, 0.20],
#     0.28: [0.38, 0.57],
#     0.30: [0.74, 0.20],
#     0.32: [0.88, 0.03],
#     0.35: [0.98, -0.10],
#     0.37: [1.03, -0.16],
#     0.38: [0.92, -0.05],
#     0.40: [0.91, -0.05],
#     0.45: [0.99, -0.17],
#     0.50: [0.97, -0.22],
#     0.55: [0.93, -0.23],
#     0.60: [0.88, -0.23],
#     0.80: [0.73, -0.22],
#     1.00: [0.64, -0.20],
#     1.5: [0.57, -0.21]
# }
#
# u2_coeff = {
#     0.20: 0.12,
#     0.22: -1.3,
#     0.245: -0.1,
#     0.265: -0.1,
#     0.28: 0.38,
#     0.30: 0.74,
#     0.32: 0.88,
#     0.35: 0.98,
#     0.37: 1.03,
#     0.38: 0.92,
#     0.40: 0.91,
#     0.45: 0.99,
#     0.50: 0.97,
#     0.55: 0.93,
#     0.60: 0.88,
#     0.80: 0.73,
#     1.00: 0.64,
# }
#
# v2_coeff = {
#     0.20: 0.33,
#     0.22: 1.6,
#     0.245: 0.85,
#     0.265: 0.20,
#     0.28: 0.57,
#     0.30: 0.20,
#     0.32: 0.03,
#     0.35: -0.10,
#     0.37: -0.16,
#     0.38: -0.05,
#     0.40: -0.05,
#     0.45: -0.17,
#     0.50: -0.22,
#     0.55: -0.23,
#     0.60: -0.23,
#     0.80: -0.22,
#     1.00: -0.2,
# }
#
#
# def interp_coeff_vals():
#     """
#     function to interpolate limb-brightening coefficient values
#
#     Returns
#     -------
#     u2_interp: interpolated u2 coefficient values
#     v2_interp: interpolated v2 coefficient values
#
#     """
#     # get coefficient lists
#     u2_coeff = get_dict('u2_coeff')
#     v2_coeff = get_dict('v2_coeff')
#
#     # get lists of mu, u2, and v2 values
#     mu_list = list(u2_coeff.keys())
#     u2_values = list(u2_coeff.values())
#     v2_values = list(v2_coeff.values())
#     mu_values = np.arange(0, 1, 0.01).tolist()
#
#     # interpolate values from dictionary
#     u2_interp1d = interp1d(mu_list, u2_values, fill_value='extrapolate')
#     u2_interp = u2_interp1d(mu_values)
#     v2_interp1d = interp1d(mu_list, v2_values, fill_value='extrapolate')
#     v2_interp = v2_interp1d(mu_values)
#
#     return u2_interp, v2_interp
#
#
# def get_limb_coeff(mu, u2_interp, v2_interp):
#
#     """
#     function to get limb darkening coefficients (Allen 1973) based on mu value
#
#     Parameters
#     ----------
#     mu: array of mu values (cos theta)
#     u2_interp: interpolated u2 coefficient values
#     v2_interp: interpolated v2 coefficient values
#
#     Returns
#     -------
#     u2: limb darkening coefficient array u2
#     v2: limb darkening coefficient array  v2
#
#     """
#
#     # get valid indices
#     use_indices = np.logical_and(mu > 0, mu != np.nan)
#
#     # get u2 values
#     u2 = np.zeros(mu.shape)
#     u2[use_indices] = [u2_interp[np.absolute(mu - mu).argmin()] for mu in mu[use_indices]]
#
#     # get v2 values
#     v2 = np.zeros(mu.shape)
#     v2[use_indices] = [v2_interp[np.absolute(mu - mu).argmin()] for mu in mu[use_indices]]
#
#     return u2, v2

