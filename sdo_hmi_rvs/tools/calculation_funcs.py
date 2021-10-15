"""
Tamar Ervin
Date: June 16, 2021

functions for calculating Solar velocity corrections
and components for derivation of SDO/HMI RVs

"""

import datetime
import numpy as np

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.coordinates import frames

from skimage.measure import label, regionprops


def map_sequence(dates_list, time_range=datetime.timedelta(seconds=6), instrument=a.Instrument.aia,
                 wavelength=a.Wavelength(171 * u.angstrom)):
    """
    function to query sunpy for images within certain time range of dates in dates list
    user specified instrument and wavelength, otherwise default values: AIA 171A

    Parameters
    ----------
    dates_list: list of dates, either datetime or strings
    time_range: plus/minus time range to search for images in comparison to desired timestamp
    instrument: Sunpy instrument of choice (AIA, HMI, LASCO, EIT)
    wavelength: desired wavelength of choice instrument

    Returns
    -------
    maps: Sunpy map sequence object

    """

    if isinstance(dates_list[0][0], str):
        datetimes = [datetime.datetime.strptime(date[0], '%Y-%m-%dT%H:%M:%S.%f') for date in dates_list]
    else:
        datetimes = dates_list

    downloaded_files = []
    for ind, datetime_object in enumerate(datetimes):
        # pull image within specified time range
        result = Fido.search(a.Time(str(datetime_object - time_range), str(datetime_object + time_range)),
                             instrument, wavelength)
        # add file to list
        downloaded_files.append(Fido.fetch(result))

    # build map sequence from files
    maps = sunpy.map.Map(downloaded_files, sequence=True)

    return maps


def rel_positions(wij, nij, rij, smap):
    """
    function to calculate pixel-wise relative positions in new coordinate frame

    Parameters
    ----------
    wij: array of westward values for image
    nij: array of northward values for image
    rij: array of radius values for image
    smap: Sunpy map object

    Returns
    -------
    deltaw: relative westward position of pixel
    deltan: relative northward position of pixel
    deltar: relative radial position of pixel
    dij: distance between pixel ij and spacecraft

    """

    # calculate relative positions of each pixel
    rsc = smap.meta['dsun_obs'] / smap.meta['rsun_ref']
    deltaw = wij
    deltan = nij
    deltar = rij - rsc
    dij = np.sqrt(deltaw ** 2 + deltan ** 2 + deltar ** 2)

    return deltaw, deltan, deltar, dij


def spacecraft_vel(deltaw, deltan, deltar, dij, vmap):
    """
    function to calculate pixel-wise spacecraft velocities for Sunpy map

    Parameters
    ----------
    deltaw: relative westward position of pixel
    deltan: relative northward position of pixel
    deltar: relative radial position of pixel
    dij: distance between pixel ij and spacecraft
    vmap: Sunpy map object (Dopplergram)

    Returns
    -------
    vsc: array of spacecraft velocities

   """

    # velocity of spacecraft relative to sun
    vscw = vmap.meta['obs_vw']
    vscn = vmap.meta['obs_vn']
    vscr = vmap.meta['obs_vr']

    # pixel-wise magnitude of spacecraft velocity
    vsc = - (deltaw * vscw + deltan * vscn + deltar * vscr) / dij

    return vsc


def solar_rot_vel(wij, nij, rij, deltaw, deltan, deltar, dij, vmap, a_parameters=[14.713, -2.396, -1.787]):
    """
    function to calculate pixel-wise velocities due to solar rotation

    Parameters
    ----------
    wij: array of westward values for image
    nij: array of northward values for image
    rij: array of radius values for image
    deltaw: relative westward position of pixel
    deltan: relative northward position of pixel
    deltar: relative radial position of pixel
    dij: distance between pixel ij and spacecraft
    vmap: Sunpy map object (Dopplergram)
    a_parameters: array of solar differential rotation parameters from Snodgrass & Ulrich (1990).

    Returns
    -------
    vrot: array of solar rotation velocities

    """

    # apply to cartesian coordinates
    x1 = wij
    y1 = nij * np.cos(np.deg2rad(vmap.meta['crlt_obs'])) + rij * np.sin(np.deg2rad(vmap.meta['crlt_obs']))
    z1 = - nij * np.sin(np.deg2rad(vmap.meta['crlt_obs'])) + rij * np.cos(np.deg2rad(vmap.meta['crlt_obs']))

    hx = x1 * np.cos(np.deg2rad(vmap.meta['crln_obs'])) + z1 * np.sin(np.deg2rad(vmap.meta['crln_obs']))
    hy = y1
    hz = -x1 * np.sin(np.deg2rad(vmap.meta['crln_obs'])) + z1 * np.cos(np.deg2rad(vmap.meta['crln_obs']))

    # apply parameters to determine vrot for given image pixel
    w = (a_parameters[0] + a_parameters[1] * ((np.sin(hy)) ** 2) + a_parameters[2] * (
            (np.sin(hy)) ** 4)) * 1. / 86400. * np.pi / 180.

    # get projection of solar rotation
    vx_rot = w * hz * vmap.meta['rsun_ref']
    vy_rot = 0.
    vz_rot = -w * hx * vmap.meta['rsun_ref']

    v1 = np.cos(np.deg2rad(vmap.meta['crln_obs'])) * vx_rot - np.sin(np.deg2rad(vmap.meta['crln_obs'])) * vz_rot
    v2 = vy_rot
    v3 = np.sin(np.deg2rad(vmap.meta['crln_obs'])) * vx_rot + np.cos(np.deg2rad(vmap.meta['crln_obs'])) * vz_rot

    # project into correct direction
    vrotw = v1
    vrotn = v2 * np.cos(np.deg2rad(vmap.meta['crlt_obs'])) - v3 * np.sin(np.deg2rad(vmap.meta['crlt_obs']))
    vrotr = v2 * np.sin(np.deg2rad(vmap.meta['crlt_obs'])) + v3 * np.cos(np.deg2rad(vmap.meta['crlt_obs']))

    # get full rotational velocity
    vrot = (deltaw * vrotw + deltan * vrotn + deltar * vrotr) / dij

    return vrot


def corrected_map(corrected_data, smap, map_type, frame=frames.HeliographicCarrington):
    """
    function to make Sunpy map object from corrected data

    Parameters
    ----------
    corrected_data: corrected velocity data
    smap: original Sunpy map object
    map_type: map type for 'content' section of fits header (string)
    frame: new rotation frame

    Returns
    -------
    corr_map: Sunpy map object with new frame information and corrected data

    """

    # build SkyCoord instance in new frame
    coord = SkyCoord(0 * u.arcsec, 0 * u.arcsec, obstime=smap.date, observer=smap.observer_coordinate,
                     frame=frame)

    # create fits header file with data and coordinate system information
    header = sunpy.map.make_fitswcs_header(corrected_data, coord)

    # update fits header with instrument and content information
    header['content'] = map_type
    header['telescop'] = smap.meta['telescop']
    header['wavelnth'] = smap.meta['wavelnth']

    # create new Sunpy map instance with corrected data
    corr_map = sunpy.map.Map(corrected_data, header)

    return corr_map


def mag_field(mu, mmap, B_noise=8, mu_cutoff=0.3):
    """
    function to correct for unsigned magnetic field strength and magnetic noise

    Parameters
    ----------
    mu: array of mu (cosine theta) values
    mmap: Sunpy map object (Magnetogram)
    B_noise: magnetic noise level in Gauss
    mu_cutoff: minimum mu cutoff value

    Returns
    -------
    Bobs: corrected observed magnetic field strength
    Br: array of corrected unsigned magnetic field strength

    """

    # get valid indices
    use_indices = np.logical_and(mu > mu_cutoff, mu != np.nan)
    mag_indices = np.logical_and(use_indices, np.abs(mmap.data) < B_noise)

    # calculate full magnetic field strength
    Bobs = mmap.data
    Br = np.full(shape=mmap.data.shape, fill_value=np.nan)
    Br[use_indices] = Bobs[use_indices] / mu[use_indices]
    Bobs[mag_indices] = 0
    Br[mag_indices] = 0

    return Bobs, Br


def mag_thresh(mu, mmap, Br_cutoff=24, mu_cutoff=0.3):
    """
    function to calculate magnetic threshold and differentiate between magnetically active regions and quiet Sun

    Parameters
    ----------
    mu: array of mu (cosine theta) values
    mmap: corrected (unsigned magnetic field) Sunpy map object (Magnetogram)
    Br_cutoff: minimum cutoff value (in Gauss) for thresholding active regions
    mu_cutoff: minimum mu cutoff value for data to ignore

    Returns
    -------
    active: weights array where active pixels are 1
    quiet: weights array where active pixels are 0

    """

    # get active region indices
    active_inds = np.where(np.abs(mmap.data) * mu > Br_cutoff)
    bad_mu = np.where(mu <= mu_cutoff)

    # make active region array
    active = np.zeros(mu.shape)
    active[active_inds] = 1.
    active[bad_mu] = 0.

    # find isolated pixels
    # get area
    y_labeled = label(active, connectivity=2, background=0)
    y_area = [props.area for props in regionprops(y_labeled)]

    # area constraint
    good_area = np.where(np.array(y_area) > 5)
    good_area = good_area[0] + 1
    active_indices = np.isin(y_labeled, good_area)

    # create weights array
    active[~active_indices] = 0

    # get quiet indices
    quiet = 1 - active

    return active, quiet


def int_thresh(map_int_cor, active, quiet):
    """
    function to do intensity thresholding and differentiate between faculae (bright) and sunspots (dark)

    Parameters
    ----------
    map_int_cor: corrected (limb-darkening) Sunpy map object (Intensitygram)
    active: weights array where active pixels are 1
    quiet: weights array where active pixels are 0

    Returns
    -------
    fac_inds: array of indices where faculae are detected
    spot_inds: array of indices where sunspots are detected

    """
    # flattened intensity data
    Iflat = map_int_cor.data

    # calculate quiet sun intensity
    int_quiet = np.nansum(Iflat * quiet) / np.nansum(quiet)

    # intensity threshold
    int_cutoff = 0.89 * int_quiet

    # get faculae
    fac_inds = np.logical_and((Iflat > int_cutoff), (active > 0.5))

    # get sunspots
    spot_inds = np.logical_and((Iflat <= int_cutoff), (active > 0.5))

    return fac_inds, spot_inds


def thresh_map(fac_inds, spot_inds):
    """
    function that creates thresholded map of sunspots (1) and faculae (2)

    Parameters
    ----------
    fac_inds: array of indices where faculae are detected
    spot_inds: array of indices where sunspots are detected

    Returns
    -------
    thresh_arr: array of values denoting faculae (1) and sunspots (2)

    """

    thresh_arr = np.full(shape=fac_inds.shape, fill_value=np.nan)
    thresh_arr[fac_inds] = 1
    thresh_arr[spot_inds] = 2

    return thresh_arr


def v_quiet(map_vel_cor, imap, quiet):
    """
    function to calculate velocity due to convective motion of quiet-Sun

    Parameters
    ----------
    map_vel_cor: corrected (velocities) Sunpy map object (Dopplergram)
    imap: UNCORRECTED Sunpy map object (Intensitygram)
    quiet: weights array where active pixels have weight = 0

    Returns
    -------
    v_quiet: quiet-Sun velocity (float)

    """

    v_quiet = np.nansum(map_vel_cor.data * imap.data * quiet) / np.nansum(
        imap.data * quiet)

    return v_quiet


def v_phot(quiet, active, Lij, vrot, imap, mu, fac_inds, spot_inds, mu_cutoff=0.3):
    """
    function to calculate photometric velocity due to rotational Doppler variation

    Parameters
    ----------
    quiet: weights array where active pixels have weight = 0
    active: weights array where active pixels have weight = 1
    Lij: limb-darkening polynomial function
    vrot: solar rotational velocity
    imap: UNCORRECTED Sunpy map object (Intensitygram)
    mu: array of mu values
    fac_inds: array of indices where faculae are detected
    spot_inds: array of indices where sunspots are detected
    mu_cutoff: minimum mu cutoff value

    Returns
    -------
    v_phot: photospheric velocity perturbation (float)

    """

    # get good mu values
    good_mu = np.where(mu > mu_cutoff)

    # calculate K scaling factor
    K = np.nansum(imap.data * Lij * quiet) / np.sum((Lij[good_mu] ** 2) * quiet[good_mu])

    # calculate photospheric velocity
    v_phot = np.nansum(np.real(vrot) * (imap.data - K * Lij) * active) / np.nansum(imap.data)

    # faculae driven photospheric velocity
    vphot_bright = np.nansum(np.real(vrot) * (imap.data - K * Lij) * fac_inds) / np.nansum(imap.data)

    # sunspots driven photospheric velocity
    vphot_spot = np.nansum(np.real(vrot) * (imap.data - K * Lij) * spot_inds) / np.nansum(imap.data)

    return v_phot, vphot_bright, vphot_spot


def v_disc(map_vel_cor, imap):
    """
    function to calculate disc-averaged velocity of Sun

    Parameters
    ----------
    map_vel_cor: corrected (velocities) Sunpy map object (Dopplergram)
    imap: UNCORRECTED Sunpy map object (Intensitygram)

    Returns
    -------
    v_disc: disc averaged velocity of Sun (float)

    """

    v_disc = np.nansum(map_vel_cor.data * imap.data) / np.nansum(imap.data)

    return v_disc


def filling_factor(mu, mmap, active, fac_inds, spot_inds, mu_cutoff=0.3):
    """
    function to calculate filling factors for faculae, sunspots, and
    total magnetically active regions
    - percentage of magnetically active pixels on the solar surface at any one time

    Parameters
    ----------
    mu: array of mu (cosine theta) values
    mmap: UNCORRECTED Sunpy map object (Magnetogram)
    active: weights array where active pixels have weight = 1
    fac_inds: array of indices where faculae are detected
    spot_inds: array of indices where sunspots are detected
    mu_cutoff: minimum mu cutoff value

    Returns
    -------
    f_bright: filling factor (%) for bright areas (faculae)
    f_spot: filling factor (%) for dark areas (sunspots)
    f_total: filling factor (%) for timestamp

    """

    # get good mu values
    good_mu = np.where(mu > mu_cutoff)

    # get number of pixels
    npix = len(mmap.data[good_mu])

    # faculae
    faculae = np.zeros(mmap.data.shape)
    faculae[fac_inds] = 1.
    f_bright = np.sum(faculae) / npix * 100

    # sunspots
    spots = np.zeros(mmap.data.shape)
    spots[spot_inds] = 1.
    f_spot = np.sum(spots) / npix * 100

    # get filling factor
    f_total = np.sum(active) / npix * 100

    return f_bright, f_spot, f_total


def unsigned_flux(map_mag_obs, imap):
    """
    calculate unsigned magnetic flux

    Parameters
    ----------
    map_mag_obs: corrected observed magnetic field strength Sunpy map object (Magnetogram)
    imap: UNCORRECTED Sunpy map object (Intensitygram)

    Returns
    -------
    unsign_flux: unsigned magnetic flux

    """

    # get data arrays
    i_data = imap.data
    m_data = map_mag_obs.data
    mabs_data = np.abs(m_data)

    # unsigned flux
    unsign_flux = np.nansum(mabs_data * i_data) / np.nansum(i_data)

    return np.abs(unsign_flux)


def area_calc(active, pixA_hem):
    """
    calculate area of active pixels for a thresholded map

    Parameters
    ----------
    active: weights array where active pixels have weight = 1
    pixA_hem: pixel areas in uHem

    Returns
    -------
    area: area of each active region weighted by its intensity

    """

    # get labeling of image
    labeled = label(active)

    # get area of active regions
    area = np.zeros(active.shape)
    props = regionprops(labeled)
    info = regionprops(labeled, pixA_hem)

    # add area to array
    for k in range(1, len(info)):
        area[props[k].coords[:, 0], props[k].coords[:, 1]] = info[k].area * info[k].mean_intensity

    return area


def area_filling_factor(active, area, mu, mmap, fac_inds, athresh=20, mu_cutoff=0.3):
    """
    calculate filling factor for regions thresholded by area
    - differentiate between large and small regions
    - differentiate between plage (large) and network (small) bright regions

    Parameters
    ----------
    active: weights array where active pixels have weight = 1
    area: area of each active region weighted by its intensity
    mu: array of mu (cosine theta) values
    mmap: UNCORRECTED Sunpy map object (Magnetogram)
    fac_inds: array of indices where faculae are detected
    athresh: area threshold value between large and small regions (in uHem)
    mu_cutoff: minimum mu cutoff value for usable data

    Returns
    -------
    f_small: filling factor (%) for small magnetically active regions
    f_large: filling factor (%) for large magnetically active regions
    f_network: filling factor (%) for network (small, bright magnetically active) regions
    f_plage: filling factor (%) for plage (large, bright magnetically active) regions
    f_nonconv: filling factor (%) for regions that do not suppress convective blueshift

    """

    # get good mu values
    good_mu = np.where(mu > mu_cutoff)

    # get number of pixels
    npix = len(mmap.data[good_mu])

    # get quiet pixels
    quiet = 1 - active

    # get filling factor for 'small' magnetic features
    small = np.zeros(mmap.data.shape)
    small_inds = np.logical_and(active > 0.5, area < athresh)
    small[small_inds] = 1.
    f_small = np.nansum(small) / npix * 100

    # get filling factor for 'large' magnetic features
    large = np.zeros(mmap.data.shape)
    large_inds = np.logical_and(active > 0.5, area > athresh)
    large[large_inds] = 1.
    f_large = np.nansum(large) / npix * 100

    # get filling factor for network (small, faculae regions)
    network = np.zeros(mmap.data.shape)
    network_inds = np.logical_and(small > 0.5, fac_inds > 0.5)
    network[network_inds] = 1.
    f_network = np.nansum(network) / npix * 100

    # get filling factor for plage (large, faculae regions)
    plage = np.zeros(mmap.data.shape)
    plage_inds = np.logical_and(large > 0.5, fac_inds > 0.5)
    plage[plage_inds] = 1.
    f_plage = np.nansum(plage) / npix * 100

    # get filling factor for small, non-convective regions
    nonconv = np.zeros(mmap.data.shape)
    nonconv_inds = np.logical_and(quiet > 0.5, small > 0.5)
    nonconv[nonconv_inds] = 1.
    f_nonconv = np.nansum(nonconv) / npix * 100

    return f_small, f_large, f_network, f_plage, f_nonconv


def area_unsigned_flux(map_mag_obs, imap, area, active, athresh=20):
    """
    calculate the magnetic flux for different regions based on area cut
    and magnetic activitiy

    Parameters
    ----------
    map_mag_obs: corrected observed magnetic field strength Sunpy map object (Magnetogram)
    imap: UNCORRECTED Sunpy map object (Intensitygram)
    area: area of each active region weighted by its intensity
    active: weights array where active pixels have weight = 1
    athresh: area threshold value between large and small regions (in uHem)

    Returns
    -------
    quiet_flux: magnetic flux of quiet-Sun regions
    ar_flux: magnetic flux of active regions
    conv_flux: magnetic flux of regions that suppress convective blueshift
    pol_flux: magnetic flux of polarized regions
    pol_conv_flux: magnetic flux of polarized regions that suppress the convective blueshift

    """

    # get data arrays
    i_data = imap.data
    m_data = map_mag_obs.data
    mabs_data = np.abs(m_data)
    quiet = 1 - active

    # get large regions array
    large = np.zeros(m_data.shape)
    large_inds = np.logical_and(active > 0.5, area > athresh)
    large[large_inds] = 1.

    # calculate relevant fluxes
    quiet_flux = np.nansum(mabs_data * i_data * quiet) / np.nansum(i_data * quiet)
    ar_flux = np.nansum(mabs_data * i_data * active) / np.nansum(i_data * active)
    conv_flux = np.nansum(mabs_data * i_data * large) / np.nansum(i_data * large)
    pol_flux = np.nansum(m_data * i_data) / np.nansum(i_data)
    pol_conv_flux = np.nansum(m_data * i_data * large) / np.nansum(i_data * large)

    return quiet_flux, ar_flux, conv_flux, pol_flux, pol_conv_flux


def area_vconv(map_vel_cor, imap, active, area, athresh=20):
    """
    calculate convective velocities for different area thresholds

    Parameters
    ----------
    map_vel_cor: corrected (velocities) Sunpy map object (Dopplergram)
    imap: UNCORRECTED Sunpy map object (Intensitygram)
    active: weights array where active pixels have weight = 1
    area: area of each active region weighted by its intensity
    athresh: area threshold value between large and small regions (in uHem)

    Returns
    -------
    vconv_quiet: convective velocity due to quiet-Sun regions
    vconv_large: convective velocity due to large active regions
    vconv_small: convective velocity due to small active regions

    """

    # get data arrays
    v_data = map_vel_cor.data
    i_data = imap.data

    # get large regions array
    large = np.zeros(v_data.shape)
    large_inds = np.logical_and(active > 0.5, area > athresh)
    large[large_inds] = 1.

    # get small regions array
    small = np.zeros(v_data.shape)
    small_inds = np.logical_and(active > 0.5, area < athresh)
    small[small_inds] = 1.

    # label the regions
    labeled = label(large)
    v_props = regionprops(labeled, v_data)
    i_props = regionprops(labeled, i_data)

    # labeled for small regions
    labeled = label(small)
    v_small = regionprops(labeled, v_data)
    i_small = regionprops(labeled, i_data)

    # get quiet regions array
    quiet = 1 - active

    # array of non-convective regions
    nonconv = np.zeros(v_data.shape)
    nonconv_inds = np.logical_or(quiet > 0.5, small > 0.5)
    nonconv[nonconv_inds] = 1.

    # velocities of non convective regions
    vel_quiet = np.nansum(v_data * i_data * quiet) / np.sum(i_data * quiet)
    vel_nonconv = np.nansum(v_data * i_data * nonconv) / np.sum(i_data * nonconv)

    # velocities of convective regions
    vconv_large = np.zeros(len(v_props))
    vconv_small = np.zeros(len(v_props))
    vconv_quiet = np.zeros(len(v_props))

    for k in range(len(v_props)):
        vconv_large[k] = v_props[k].area * (v_props[k].mean_intensity - vel_quiet) * i_props[k].mean_intensity
        vconv_small[k] = v_small[k].area * (v_small[k].mean_intensity - vel_quiet) * i_small[k].mean_intensity
        vconv_quiet[k] = v_props[k].area * (v_props[k].mean_intensity - vel_nonconv) * i_props[k].mean_intensity

    # convective velocity of quiet regions
    vconv_quiet = np.nansum(vconv_quiet) / np.sum(i_data)

    # convective velocity of large regions
    vconv_large = np.nansum(vconv_large) / np.sum(i_data)

    # convective velocity of small regions
    vconv_small = np.nansum(vconv_small) / np.sum(i_data)

    return vconv_quiet, vconv_large, vconv_small



