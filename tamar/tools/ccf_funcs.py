"""
Tamar Ervin
Date: July 22, 2021

functions for calculating
relevant ccf values

"""
from __future__ import print_function, division

import numpy as np
from scipy import interpolate

from scipy import constants
from math import ceil
from astropy.modeling import models, fitting
import logging

from astropy.coordinates import get_body_barycentric_posvel
import astropy.units as u
from barycorrpy import PINT_erfautils as PINT
from barycorrpy.utils import CalculatePositionVector
from barycorrpy.PhysicalConstants import *
from scipy.spatial.transform import Rotation as R
from scipy.integrate import quad


def shift_ccf(vel_arr, ccf, rv_shift):
    """
    function to shift CCF such that it is centered around
    a velocity of 0 km/s
    - shift is based on radial velocity
        - can either be SDO derived or ground based

    Parameters
    ----------
    vel_arr: original velocity array
    cff: ccf that will be shifted
    rv_shift

    Returns
    -------
    ccf_shift: shifted ccf

    """

    # calculate array of velocities to shift by
    shift_arr = np.linspace(vel_arr[0], vel_arr[-1], len(ccf))

    # account for shift from 0 due to radial velocity
    shift_vel = shift_arr - rv_shift

    # interpolate to get new ccf values
    f = interpolate.interp1d(shift_vel, ccf, kind='cubic', fill_value='extrapolate')

    # apply interpolation function to velocity array
    ccf_shift = f(shift_arr)

    return ccf_shift


def CalculateFWHMDifference_SolarRotation_Ecliptic(loc, JDUTC):
    """
    Calculate the difference between the Observed Solar FWHM and Sidereal Solar FWHM
    Based on Collier Cameron et al. (2019)

    Parameters
    ----------
    loc: Astropy Earth Location object. https://docs.astropy.org/en/stable/api/astropy.coordinates.EarthLocation.html
    JDUTC: Astropy Time Object

    Returns
    -------
    Delta:  F_obs**2 - F_sid**2 [(km/s)^2]

    """

    # parameters
    ephemeris = 'de430'

    # Convert times to obtain TDB and TT
    JDTDB = JDUTC.tdb
    JDTT = JDUTC.tt

    ################################
    ###### EARTH EPHEMERIS ############
    ################################
    ##### NUTATION, PRECESSION, ETC. #####
    # Observatory position wrt Geocenter
    r_pint, v_pint = PINT.gcrs_posvel_from_itrf(loc, JDUTC, JDTT)
    r_eci = r_pint[0]  # [m]
    v_eci = v_pint[0]  # [m/s]
    ##### EPHEMERIDES #####
    earth_geo = get_body_barycentric_posvel('earth', JDTDB, ephemeris=ephemeris)  # [km]
    r_geo = np.reshape(earth_geo[0].xyz.value * 1000., 3)  # [m]
    v_geo = np.reshape(earth_geo[1].xyz.value * 1000. / 86400., 3)  # [m/s]
    PosVector_EarthSSB = r_eci + r_geo  # [m]
    # Relativistic Addition of Velocities
    VelVector_EarthSSB = (v_eci + v_geo) / (1. + np.sum(v_eci * v_geo) / c ** 2)  # [m/s]

    ################################
    ###### SOLAR EPHEMERIS ############
    ################################
    solar_ephem = get_body_barycentric_posvel('sun', JDTDB, ephemeris=ephemeris)
    PosVector_SolSSB = np.reshape(solar_ephem[0].xyz.value * 1000., 3)  # [m]
    VelVector_SolSSB = np.reshape(solar_ephem[1].xyz.value * 1000. / 86400., 3)  # [m/s]

    ################################
    ####EQUATORIAL COORD VECTORS ####
    ################################
    PosVector_EarthSol, PosMag_EarthSol, PosHat_EarthSol = CalculatePositionVector(r1=PosVector_EarthSSB,
                                                                                   r2=PosVector_SolSSB)
    VelVector_EarthSol = (VelVector_EarthSSB - VelVector_SolSSB) / (
                1. + np.sum(VelVector_SolSSB * VelVector_EarthSSB) / c ** 2)
    OmegaVector = np.cross(PosVector_EarthSol, VelVector_EarthSol) / (PosMag_EarthSol ** 2)
    ################################
    ################################
    DeltaCentury = (JDUTC.datetime.year - 2000) / 100
    # https://www2.mps.mpg.de/homes/fraenz/systems/systems3art.pdf
    # Eqn 14
    SolarInclination = 7.25
    SolarLongitude = 75.76 + 1.397 * DeltaCentury
    EclipticEpsilon = 23.44
    # Need to perform Extrinsic Euler Rotations to rotate from equatorial to ecliptic
    REcliptic = R.from_euler("X", -EclipticEpsilon, degrees=True).as_matrix()
    # Intrinsic rotation to go from solar axis to ecliptic
    RObliquity = R.from_euler("xz", [SolarInclination, SolarLongitude], degrees=True).as_matrix()

    ################################
    ####ROTATED COORD VECTORS ####
    ################################
    RotatedPositionVector = np.matmul(REcliptic, PosVector_EarthSol)
    RotatedPositionHat = RotatedPositionVector / np.linalg.norm(RotatedPositionVector)
    OmegaEarthVectorRotated = np.matmul(REcliptic, OmegaVector)
    ################################
    ################################
    OmegaSolVector = np.array([0, 0, 2.972 * 1e-6])
    OmegaSolHat = np.array([0, 0, 1])
    # Rotated Solar rotation vector to ecliptic plane
    OmegaSolVector = np.matmul(RObliquity, OmegaSolVector)
    OmegaSolHat = OmegaSolVector / np.linalg.norm(OmegaSolVector)
    sini = np.sqrt(1 - np.matmul(OmegaSolHat, RotatedPositionHat) ** 2)
    Gamma = 1.04
    DeltaOmega = OmegaSolVector - OmegaEarthVectorRotated
    Delta = ((Gamma * ac.R_sun.to(u.km).value) ** 2) * (
                np.matmul(DeltaOmega, DeltaOmega) * sini * sini - np.matmul(OmegaSolVector, OmegaSolVector))

    return Delta


def normalize_residual(residual_ccf):
    """
    function to normalize the residual CCFs
    Parameters
    ----------
    residual_ccf: residual CCF where the quiet-Sun template CCF is subtracted from the active day CCF

    Returns
    -------
    ccf: normalized residual ccf

    """

    ccf = residual_ccf - np.median(residual_ccf)
    ccf /= np.std(residual_ccf)

    return ccf


def fit_gaussian_to_ccf(velocity_loop, ccf, rv_guess, velocity_halfrange_to_fit=100.0):
    """
    Fit a Gaussian to the cross-correlation function: Arpita Roy
        - Note that median value of the CCF is subtracted before
    the fit, since astropy has trouble otherwise.
        - Analyses of CCF absolute levels must be performed separately.
    Parameters
    ----------
    velocity_loop: array of velocity values corresponding with CCF
    ccf: cross correlation function
    rv_guess: estimate of RV (from NEID?)
    velocity_halfrange_to_fit: range around which we fit the CCF

    Returns
    -------
    gaussian_fit: fit of CCF
    g_x: velocity loop for range around CCF peak
    g_y: CCF for range around CCF peak
    rv_mean: mean of fitted Gaussian - calculated RV value

    """

    FIT_G = fitting.LevMarLSQFitter()
    g_init = models.Gaussian1D(amplitude=-1e7, mean=rv_guess, stddev=2.5)

    ccf = ccf - np.nanmedian(ccf)

    # First look for CCF peak around user-supplied stellar gamma velocity
    gaussian_fit_init = FIT_G(g_init, velocity_loop, ccf)

    # If the peak doesn't look right
    if gaussian_fit_init.amplitude > 0 \
            or gaussian_fit_init.amplitude > (0 - np.std(ccf)) \
            or gaussian_fit_init.stddev < 2 \
            or gaussian_fit_init.mean > np.max(velocity_loop) \
            or gaussian_fit_init.mean < np.min(velocity_loop):

        # Try fitting minimum of CCF
        crude_guess = velocity_loop[np.argmin(ccf)]
        g_init = models.Gaussian1D(amplitude=-1e7, mean=crude_guess, stddev=2.5)
        gaussian_fit_init = FIT_G(g_init, velocity_loop, ccf)

        # If it still doesn't conform to expectation
        if gaussian_fit_init.amplitude > 0 \
                or gaussian_fit_init.amplitude > (0 - np.std(ccf)) \
                or gaussian_fit_init.stddev < 1 \
                or gaussian_fit_init.mean > np.max(velocity_loop) \
                or gaussian_fit_init.mean < np.min(velocity_loop):
            logging.error('No significant peak found in CCF')

    g_init = models.Gaussian1D(
        amplitude=gaussian_fit_init.amplitude,
        mean=gaussian_fit_init.mean,
        stddev=gaussian_fit_init.stddev)
    rv_guess = gaussian_fit_init.mean

    # Fit smaller range around CCF peak
    i_fit = ((velocity_loop >= rv_guess - velocity_halfrange_to_fit) &
             (velocity_loop <= rv_guess + velocity_halfrange_to_fit))
    g_x = velocity_loop[i_fit]
    g_y = ccf[i_fit]
    gaussian_fit = FIT_G(g_init, g_x, g_y)
    rv_mean = gaussian_fit.mean.value

    return gaussian_fit, g_x, g_y, rv_mean


def gauss(x, gaussian_fit):
    """
    function to define gaussian function to integrate area

    Parameters
    ----------
    x: x array to integrate over
    gaussian_fit: gaussian fitted ccf (Gaussian1D object)

    Returns
    -------
    g: gaussian function

    """

    # gaussian function built off fit parameters from astropy's Guassian 1D
    g = gaussian_fit.amplitude * np.exp(- ((x - gaussian_fit.mean) ** 2) / (2 * (gaussian_fit.stddev ** 2)))

    # invert this if needed such that we get a positive area
    if np.min(g) < 0:
        g = - g

    return g


def integrated_area(gaussian_fit):
    """
    function to calculate integrated area of fitting gaussian ccf

    Parameters
    ----------
    gaussian_fit: gaussian fitted ccf (Gaussian1D object)

    Returns
    -------
    int_area: integrated area of fitting ccf
    int_error: error from integrated area of ccf

    """

    # integrate to get the area under the curve!! yay exciting
    int_area, int_error = quad(gauss(gaussian_fit), gaussian_fit.bounding_box[0], gaussian_fit.bounding_box[1])

    return int_area, int_error


def bisector():
    # get left point of choice
    # np.argwhere( )
    left = 0

    # get corresponding right point... i guess draw horizontal line??
    # probably some function to do this
    right = 1

    # get midpoint -- lolz average
    mdpt = np.average(left, right)

    # also at somepoint calculate the velocity span...this is in km/s
    # v_span = v_top - v_bottom

    # return v_top, v_bottom, v_span
    return None