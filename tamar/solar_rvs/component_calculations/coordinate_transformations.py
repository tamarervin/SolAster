"""
Tamar Ervin
Date: July 8, 2021

updates to coordinate transformations to work with
solar rotational velocity calculation
"""

import numpy as np
import datetime
import astropy.units as u
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


def get_scales_from_map(smap):
    """
    compute the solar X and solar Y 1D pixel scale arrays from a sunpy map object
    - If the image has been rotated to solar north up, then the x and y scales will
      be in the helioprojective cartesian system with units of [Rs].

    Parameters
    ----------
    smap: Sunpy map object

    Returns
    -------
    x: array of x coordinates in helioprojective cartesian system
    y: array of y coordinates in helioprojective cartesian system

    """

    # get the x, y, pixel coordinates in 1D by going along the diagonal (assume its square)
    npix = int(smap.dimensions[0].value)
    inds = np.arange(npix) * u.pix
    xylocs = smap.pixel_to_world(inds, inds)

    # Tx and Ty are the arcseconds as astropy quantities
    x_rs = xylocs.Tx / smap.rsun_obs
    y_rs = xylocs.Ty / smap.rsun_obs

    # convert them to floating point
    x = np.float32(x_rs.value)
    y = np.float32(y_rs.value)

    return x, y


def map_to_image_rot_mat(obsv_lon, obsv_lat):
    """
    calculate rotation matrix for coordinate transformation

    Parameters
    ----------
    obsv_lon: longitude of observing telescope
    obsv_lat: latitude of observing telescope

    Returns
    -------
    tot_rot: rotation matrix

    """

    # rotate phi (about map z-axis. observer phi goes to -y)
    del_phi = -obsv_lon * np.pi / 180 - np.pi / 2
    rot1 = np.array([[np.cos(del_phi), -np.sin(del_phi), 0.],
                     [np.sin(del_phi), np.cos(del_phi), 0.], [0., 0., 1.], ])

    # rotate theta (about x-axis. observer theta goes to +z)
    del_theta = obsv_lat * np.pi / 180 - np.pi / 2
    rot2 = np.array([[1., 0., 0.], [0., np.cos(del_theta), -np.sin(del_theta)],
                     [0., np.sin(del_theta), np.cos(del_theta)]])

    tot_rot = np.matmul(rot2, rot1)

    return tot_rot


def image_grid_to_cr(x, y, R0=1.01, obsv_lat=0, obsv_lon=0, outside_map_val=np.nan):
    """
    function to transform vector coordinate pairs in solar radii units and the observer angles to map coordinates

    Parameters
    ----------
    x: array of helioprojective cartesian x values
    y: array of helioprojective cartesian y values
    R0: observation radius  # TODO: determine the correct R0 value to use -- rsun1 = mag_map_rot.meta['rsun_obs'] * 760e3 / mag_map_rot.meta['rsun_ref']
    obsv_lon: longitude of observing telescope
    obsv_lat: latitude of observing telescope
    outside_map_val: value for array elements outside of the solar limb

    Returns
    -------
    cr_theta_all: array of theta values in carrington coordinates
    cr_phi_all: array of phi values in carrington coordinates
    image_mu: array of mu values (cos theta)

    """

    # for images, we assume that the z-axis is perpendicular to the image plane, in the direction
    # of the observer, and located at the center of the image.

    # mask points outside of R0
    use_index = x ** 2 + y ** 2 <= R0 ** 2
    use_x = x[use_index]
    use_y = y[use_index]

    # Find z coord (we can assume it is in the positive direction)
    # to be numerically equivalent to the use_index definition, change to this:
    use_z = np.sqrt(R0 ** 2 - (use_x ** 2 + use_y ** 2))

    # Calc image_theta, image_phi, and image_mu
    image_mu = np.full(x.shape, 0)
    use_theta = np.arccos(use_z / R0)
    image_mu[use_index] = np.cos(use_theta)

    # generate map-to-image rotation matrix
    rot_mat = map_to_image_rot_mat(obsv_lon, obsv_lat)
    # invert/transpose for image-to-map rotation matrix
    rev_rot = rot_mat.transpose()
    # construct coordinate array
    coord_array = np.array([use_x, use_y, use_z])
    # apply rotation matrix to coordinates
    map3D_coord = np.matmul(rev_rot, coord_array)

    # Occasionally numeric error from the rotation causes a z magnitude to be greater than R0
    num_err_z_index = np.abs(map3D_coord[2, :]) > R0
    map3D_coord[2, num_err_z_index] = np.sign(map3D_coord[2, num_err_z_index]) * R0

    # Convert map cartesian to map theta and phi
    cr_theta = np.arccos(map3D_coord[2, :] / R0)
    cr_phi = np.arctan2(map3D_coord[1, :], map3D_coord[0, :])
    cr_r = np.sqrt(map3D_coord[0, :] ** 2 + map3D_coord[1, :] ** 2 + map3D_coord[2, :] ** 2)

    # Change phi range from [-pi,pi] to [0,2pi] -- TODO: i no longer think this should be done
    # neg_phi = cr_phi < 0
    # cr_phi[neg_phi] = cr_phi[neg_phi] + 2 * np.pi

    cr_theta_all = np.full(x.shape, outside_map_val)
    cr_phi_all = np.full(x.shape, outside_map_val)
    cr_r_all = np.full(x.shape, outside_map_val)

    cr_theta_all[use_index] = cr_theta
    cr_phi_all[use_index] = cr_phi
    cr_r_all[use_index] = cr_r

    return cr_theta_all, cr_phi_all, cr_r_all, image_mu


def get_coordinates(x, y, smap, R0=1.01, outside_map_val=np.nan):
    """
    Calculate relevant mapping information for each pixel.

    Parameters
    ----------
    x: array of x coordinates in helioprojective cartesian system
    y: array of y coordinates in helioprojective cartesian system
    smap: Sunpy map object
    R0: observation radius
    outside_map_val: value for array elements outside of the solar limb

    Returns
    -------
    lat: Carrington latitude
    lon: Carrington longitude
    mu: cosine of the center to limb angle

    """

    # create grids for x, y arrays
    x_mat, y_mat = np.meshgrid(x, y)
    x_vec = x_mat.flatten(order="C")
    y_vec = y_mat.flatten(order="C")

    # calculate theta, phi, and mu values in Carrington frame
    cr_theta_all, cr_phi_all, cr_r_all, image_mu = image_grid_to_cr(x_vec, y_vec, R0=R0, obsv_lat=smap.meta['crlt_obs'],
                                                    obsv_lon=smap.meta['crln_obs'], outside_map_val=outside_map_val)

    # reshape arrays to match map dimensions
    cr_theta = cr_theta_all.reshape(smap.data.shape, order="C")
    cr_phi = cr_phi_all.reshape(smap.data.shape, order="C")
    cr_r = cr_r_all.reshape(smap.data.shape, order="C")
    image_mu = image_mu.reshape(smap.data.shape, order="C")

    # calculate latitude and longitude from theta and phi values
    lat = cr_theta - np.pi / 2.
    lon = cr_phi
    r = cr_r
    mu = image_mu

    return lat, lon, r, mu


def coord_trans(smap, R0=1.01, outside_map_val=np.nan):
    """
    function for calculation of coordinate transformations of HMI images

    Parameters
    ----------
    smap: Sunpy map object
    R0: observation radius
    outside_map_val: value for array elements outside of the solar limb

    Returns
    -------
    theta: Carrington latitude
    phi: Carrington longitude
    mu: cosine of the center to limb angle

    """

    x, y = get_scales_from_map(smap)
    lat, lon, r, mu = get_coordinates(x, y, smap, R0, outside_map_val=outside_map_val)

    return lat, lon, r, mu


lat, lon, r, mymu = coord_trans(smap, R0=smap.meta['rsun_obs'] * 760e3 / smap.meta['rsun_ref'])

import matplotlib.pyplot as plt
image_mu = np.where(image_mu == np.nan, 0, image_mu)
plt.imshow(image_mu)
plt.title('my mu')
plt.colorbar()
plt.show()

import tamar.tools.coord_funcs as ctfuncs
px, py, p, r, d, mu = ctfuncs.coordinates(smap)
rw_obs, rn_obs, rr_obs = ctfuncs.vel_coords(px, py, p, r, smap)
plt.imshow(mu)
plt.title('mu')
plt.colorbar()
plt.show()
