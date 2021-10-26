"""
Tamar Ervin
Date: June 24, 2021

Coordinate transformation functions to transform from the Helioprojective
Cartesian to Heliographic Carrington Coordinates.

"""

import scipy
import numpy as np


def get_map_scales(smap):
    """
    compute the solar X and solar Y 2D pixel scale arrays from a sunpy map object
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

    header = smap.fits_header

    # pixel array
    axis1 = np.linspace(1, header['NAXIS1'], header['NAXIS1']) - header['CRPIX1']
    axis2 = np.linspace(1, header['NAXIS2'], header['NAXIS2']) - header['CRPIX2']

    # pixel offset from center of image, x,y
    x, y = np.meshgrid(axis1, axis2)

    return x, y


def coordinates(smap):
    """
    calculate array of mu values and cartesian
    coordinates for image

    Parameters
    ----------
    smap: Sunpy map object

    Returns
    -------
    x: array of x coordinates in helioprojective cartesian system
    y: array of y coordinates in helioprojective cartesian system
    pd: array of pixel distances
    r: array of solar radius values
    d: observer distance from sun in solar radii
    mu: array of mu values


    """

    # distance of observer to sun in solar radii
    d = smap.fits_header['DSUN_OBS'] / smap.fits_header['RSUN_REF']

    # focal length in pixels
    f = 180. * 3600. / np.pi / smap.fits_header['CDELT1']

    # get cartesian x, y map scales
    x, y = get_map_scales(smap)

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

    return x, y, pd, r, d, mu


def vel_coords(x, y, pd, r, smap):
    """
    calculate coordinate transformation to heliographic Carrington coordinates

    Parameters
    ----------
    x: array of x coordinates in helioprojective cartesian system
    y: array of y coordinates in helioprojective cartesian system
    pd: array of pixel distances
    r: array of solar radius values
    smap: Sunpy map object

    Returns
    -------
    wij: array of pixel coordinates relative to solar center in westward direction
    nij: array of pixel coordinates relative to solar center in northward direction
    rij: array of pixel coordinates relative to solar center in radial direction

    """

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

    return wij, nij, rij


def fix_mu(mu, smaps, mu_cutoff=0.3):
    """
    function to remove pixel values where mu is less than 0.1

    Parameters
    ----------
    mu: cosine of the center to limb angle
    smaps: list of Sunpy map object
    mu_cutoff: minimum mu cutoff value

    Returns
    -------
    mu: corrected cosine of the center to limb angle
    smaps: corrected Sunpy map objects

    """

    # remove pixel values where mu < mu_cut (0.3)
    bad_inds = np.where(mu <= mu_cutoff)

    # fix arrays
    for smap in smaps:
        smap.data[bad_inds] = 0

    return smaps


def pix_area_hem(wij, nij, rij, smap):
    """
    calculate the area of each pixel in uHem for use in area thresholding
    for convective velocity and identification of solar regions

    Parameters
    ----------
    wij: array of pixel coordinates relative to solar center in westward direction
    nij: array of pixel coordinates relative to solar center in northward direction
    rij: array of pixel coordinates relative to solar center in radial direction
    smap

    Returns
    -------
    pixA_hem: pixel areas in uHem

    """
    # get x and y arrays
    x, y = get_map_scales(smap)

    # apply to cartesian coordinates
    x1 = wij
    y1 = nij * np.cos(np.deg2rad(smap.meta['crlt_obs'])) + rij * np.sin(np.deg2rad(smap.meta['crlt_obs']))
    z1 = - nij * np.sin(np.deg2rad(smap.meta['crlt_obs'])) + rij * np.cos(np.deg2rad(smap.meta['crlt_obs']))

    # get solar latitude and longitude
    solar_long = (np.arctan(x1 / z1)).real
    solar_long = -abs(solar_long) * np.sign(x)
    y1 = y1 + 0.J
    solar_lat = (np.arcsin(y1)).real
    d_long = np.diff(solar_long, 1, 1)
    d_lat = np.diff(solar_lat, 1, 0)

    # interpolate differential longitude and latitude to get same area dimensions as original grid
    newxx = x[:, 0:-1:] + 1 / 2 * np.diff(x, axis=1)
    newyx = y[:, 0:-1:] + 1 / 2 * np.diff(y, axis=1)
    newxy = x[0:-1:, :] + 1 / 2 * np.diff(x, axis=0)
    newyy = y[0:-1:, :] + 1 / 2 * np.diff(y, axis=0)

    # interpolate latitude array
    interp = scipy.interpolate.interp2d(newxy[0], newyy[:, 0], d_lat)
    d_lat = interp(x[0], y[:, 0])

    # interpolate longitude array
    interp = scipy.interpolate.interp2d(newxx[0], newyx[:, 0], d_long)
    d_long = interp(x[0], y[:, 0])

    # calculate co-latitude area in uHem
    pixA_hem = np.sin(np.pi / 2 - solar_lat) * abs(d_long) * abs(d_lat) / (2 * np.pi) * 1e6

    return pixA_hem


# defunct

# def get_mu(x, y, smap):
#     """
#     calculate array of mu values for image
#
#     Parameters
#     ----------
#     x: array of x coordinates in helioprojective cartesian system
#     y: array of y coordinates in helioprojective cartesian system
#     smap: Sunpy map object
#
#     Returns
#     -------
#     mu: array of mu values
#     r: array of solar radius values
#
#     """
#     # pixelwise distance
#     pd = np.sqrt(x ** 2 + y ** 2)
#
#     # distance of observer to sun in solar radii
#     d = smap.fits_header['DSUN_OBS'] / smap.fits_header['RSUN_REF']
#
#     # focal length in pixels
#     f = 180. * 3600. / np.pi / smap.fits_header['CDELT1']
#
#     # distance (in solar radii) to pixel
#     r = f * f * pd * pd + pd ** 4 - d * d * pd ** 4 + 0.J
#     r = (d * f * pd - np.sqrt(r)) / (f * f + pd * pd)
#
#     # we only want the real part
#     r = r.real
#
#     # get mu array
#     r = 1 - r ** 2 + 0.J
#     cos_alpha = (np.sqrt(r)).real
#     sin_alpha = r.real
#     cos_theta = ((d - np.sqrt(r)) / np.sqrt(r ** 2 + (d - np.sqrt(r)) ** 2)).real
#     sin_theta = (np.sqrt(1 - cos_theta ** 2)).real
#     mu = cos_alpha * cos_theta - sin_alpha * sin_theta
#
#     return mu, r
#
#
# def coord_trans(x, y, r, smap):
#     """
#     calculate coordinate transformation to heliographic Carrington coordinates
#
#     Parameters
#     ----------
#     x: array of x coordinates in helioprojective cartesian system
#     y: array of y coordinates in helioprojective cartesian system
#     r: array of solar radius values
#     smap: Sunpy map object
#
#     Returns
#     -------
#     wij: array of pixel coordinates relative to solar center in westward direction
#     nij: array of pixel coordinates relative to solar center in northward direction
#     rij: array of pixel coordinates relative to solar center in radial direction
#
#     """
#     # pixelwise distance
#     pd = np.sqrt(x ** 2 + y ** 2)
#
#     # transform each pixel to get into Heliographic CR Coordinates
#     dw = x * np.cos(np.deg2rad(smap.fits_header['CROTA2'])) + y * np.sin(np.deg2rad(smap.fits_header['CROTA2']))
#     dn = - x * np.sin(np.deg2rad(smap.fits_header['CROTA2'])) + y * np.cos(np.deg2rad(smap.fits_header['CROTA2']))
#
#     # get cartesian coordinates for velocity calculations
#     wij = r * dw / pd
#     nij = r * dn / pd
#     rij = (np.sqrt(1 - r ** 2 + 0.J)).real
#
#     return wij, nij, rij
#
#
# def get_scales_from_map(smap):
#     """
#     compute the solar X and solar Y 1D pixel scale arrays from a sunpy map object
#     - If the image has been rotated to solar north up, then the x and y scales will
#       be in the helioprojective cartesian system with units of [Rs].
#
#     Parameters
#     ----------
#     smap: Sunpy map object
#
#     Returns
#     -------
#     x: array of x coordinates in helioprojective cartesian system
#     y: array of y coordinates in helioprojective cartesian system
#
#     """
#
#     # get the x, y, pixel coordinates in 1D by going along the diagonal (assume its square)
#     npix = int(smap.dimensions[0].value)
#     inds = np.arange(npix) * u.pix
#     xylocs = smap.pixel_to_world(inds, inds)
#
#     # Tx and Ty are the arcseconds as astropy quantities
#     x_rs = xylocs.Tx / smap.rsun_obs
#     y_rs = xylocs.Ty / smap.rsun_obs
#
#     # convert them to floating point
#     x = np.float32(x_rs.value)
#     y = np.float32(y_rs.value)
#
#     return x, y
#
#
# def map_to_image_rot_mat(obsv_lon, obsv_lat):
#     """
#     calculate rotation matrix for coordinate transformation
#
#     Parameters
#     ----------
#     obsv_lon: longitude of observing telescope
#     obsv_lat: latitude of observing telescope
#
#     Returns
#     -------
#     tot_rot: rotation matrix
#
#     """
#
#     # rotate phi (about map z-axis. observer phi goes to -y)
#     del_phi = -obsv_lon * np.pi / 180 - np.pi / 2
#     rot1 = np.array([[np.cos(del_phi), -np.sin(del_phi), 0.],
#                      [np.sin(del_phi), np.cos(del_phi), 0.], [0., 0., 1.], ])
#
#     # rotate theta (about x-axis. observer theta goes to +z)
#     del_theta = obsv_lat * np.pi / 180 - np.pi / 2
#     rot2 = np.array([[1., 0., 0.], [0., np.cos(del_theta), -np.sin(del_theta)],
#                      [0., np.sin(del_theta), np.cos(del_theta)]])
#
#     tot_rot = np.matmul(rot2, rot1)
#
#     return tot_rot
#
#
# def image_grid_to_cr(x, y, R0=1.01, obsv_lat=0, obsv_lon=0, outside_map_val=np.nan):
#     """
#     function to transform vector coordinate pairs in solar radii units and the observer angles to map coordinates
#
#     Parameters
#     ----------
#     x: array of helioprojective cartesian x values
#     y: array of helioprojective cartesian y values
#     R0: observation radius
#     obsv_lon: longitude of observing telescope
#     obsv_lat: latitude of observing telescope
#     outside_map_val: value for array elements outside of the solar limb
#
#     Returns
#     -------
#     cr_theta_all: array of theta values in carrington coordinates
#     cr_phi_all: array of phi values in carrington coordinates
#     image_mu: array of mu values (cos theta)
#
#     """
#
#     # for images, we assume that the z-axis is perpendicular to the image plane, in the direction
#     # of the observer, and located at the center of the image.
#
#     # mask points outside of R0
#     use_index = x ** 2 + y ** 2 <= R0 ** 2
#     use_x = x[use_index]
#     use_y = y[use_index]
#
#     # Find z coord (we can assume it is in the positive direction)
#     # to be numerically equivalent to the use_index definition, change to this:
#     use_z = np.sqrt(R0 ** 2 - (use_x ** 2 + use_y ** 2))
#
#     # Calc image_theta, image_phi, and image_mu
#     image_mu = np.full(x.shape, outside_map_val)
#     use_theta = np.arccos(use_z / R0)
#     image_mu[use_index] = np.cos(use_theta)
#
#     # generate map-to-image rotation matrix
#     rot_mat = map_to_image_rot_mat(obsv_lon, obsv_lat)
#     # invert/transpose for image-to-map rotation matrix
#     rev_rot = rot_mat.transpose()
#     # construct coordinate array
#     coord_array = np.array([use_x, use_y, use_z])
#     # apply rotation matrix to coordinates
#     map3D_coord = np.matmul(rev_rot, coord_array)
#
#     # Occasionally numeric error from the rotation causes a z magnitude to be greater than R0
#     num_err_z_index = np.abs(map3D_coord[2, :]) > R0
#     map3D_coord[2, num_err_z_index] = np.sign(map3D_coord[2, num_err_z_index]) * R0
#
#     # Convert map cartesian to map theta and phi
#     cr_theta = np.arccos(map3D_coord[2, :] / R0)
#     cr_phi = np.arctan2(map3D_coord[1, :], map3D_coord[0, :])
#     cr_r = np.sqrt(map3D_coord[0, :] ** 2 + map3D_coord[1, :] ** 2 + map3D_coord[2, :] ** 2)
#
#     # Change phi range from [-pi,pi] to [0,2pi]
#     # neg_phi = cr_phi < 0
#     # cr_phi[neg_phi] = cr_phi[neg_phi] + 2 * np.pi
#
#     cr_theta_all = np.full(x.shape, outside_map_val)
#     cr_phi_all = np.full(x.shape, outside_map_val)
#     cr_r_all = np.full(x.shape, outside_map_val)
#
#     cr_theta_all[use_index] = cr_theta
#     cr_phi_all[use_index] = cr_phi
#     cr_r_all[use_index] = cr_r
#
#     return cr_theta_all, cr_phi_all, cr_r_all, image_mu
#
#
# def get_coordinates(x, y, smap, R0=1.01, outside_map_val=np.nan):
#     """
#     Calculate relevant mapping information for each pixel.
#
#     Parameters
#     ----------
#     x: array of x coordinates in helioprojective cartesian system
#     y: array of y coordinates in helioprojective cartesian system
#     smap: Sunpy map object
#     R0: observation radius
#     outside_map_val: value for array elements outside of the solar limb
#
#     Returns
#     -------
#     lat: Carrington latitude
#     lon: Carrington longitude
#     mu: cosine of the center to limb angle
#
#     """
#
#     # create grids for x, y arrays
#     x_mat, y_mat = np.meshgrid(x, y)
#     x_vec = x_mat.flatten(order="C")
#     y_vec = y_mat.flatten(order="C")
#
#     # calculate theta, phi, and mu values in Carrington frame
#     cr_theta_all, cr_phi_all, cr_r_all, image_mu = image_grid_to_cr(x_vec, y_vec, R0=R0, obsv_lat=smap.meta['crlt_obs'],
#                                                     obsv_lon=smap.meta['crln_obs'], outside_map_val=outside_map_val)
#
#     # reshape arrays to match map dimensions
#     cr_theta = cr_theta_all.reshape(smap.data.shape, order="C")
#     cr_phi = cr_phi_all.reshape(smap.data.shape, order="C")
#     cr_r = cr_r_all.reshape(smap.data.shape, order="C")
#     image_mu = image_mu.reshape(smap.data.shape, order="C")
#
#     # calculate latitude and longitude from theta and phi values
#     lat = cr_theta - np.pi / 2.
#     lon = cr_phi
#     r = cr_r
#     mu = image_mu
#
#     return lat, lon, r, mu







