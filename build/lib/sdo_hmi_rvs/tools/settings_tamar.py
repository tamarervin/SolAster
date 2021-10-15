"""
Tamar Ervin
Date: June 14, 2021

Settings file for paths
Update to work with your system.
"""

import os


class BaseDir:
    """
    Base directory. Update this and follow same file structure and all will be well :)
    """
    BASE_DIR = '/Volumes/Portable'


class CsvDir:
    """
    CSV directories
    """
    CSV_DIR = os.path.join(BaseDir.BASE_DIR, 'csv_files')
    MILBOURNE = os.path.join(CSV_DIR, 'milbourne')
    NEID_DATA = os.path.join(CSV_DIR, 'neid', 'data')
    NEID_CALC = os.path.join(CSV_DIR, 'neid', 'calcs')
    NEID_BAD_DATES = os.path.join(CSV_DIR, 'neid', 'bad_dates')
    VENUS = os.path.join(CSV_DIR, 'venus')
    MERCURY = os.path.join(CSV_DIR, 'mercury')
    ACTIVE = os.path.join(CSV_DIR, 'active')
    CCFS = os.path.join(CSV_DIR, 'ccfs')
    NEID_SOLAR = os.path.join(BaseDir.BASE_DIR, 'mnt', 'grinnell', 'NEID', 'NEIDdata', 'DRPProcessed', 'DRPmaster', 'solar')
    WEIGHTS = os.path.join(BaseDir.BASE_DIR, 'NEID_Solar_analysis', 'example_notebooks_and_data', 'ord_weights_test.npy')
    FSR_MASK = os.path.join(BaseDir.BASE_DIR, 'NEID_Solar_analysis', 'NEIDcode', 'masks', 'neidMaster_FSR_Mask20210331_v001.fits')
    NEID_HOUR = os.path.join(BaseDir.BASE_DIR, 'neid_data')
    LEVEL0 = os.path.join(BaseDir.BASE_DIR, 'level0')


class ImgDir:
    """
    Image directories
    """
    IMG_DIR = os.path.join(BaseDir.BASE_DIR, 'images')
    NEID_IMG = os.path.join(IMG_DIR, 'neid')
    MERCURY_IMG = os.path.join(IMG_DIR, 'mercury')
    VENUS_IMG = os.path.join(IMG_DIR, 'venus')
    CCF_IMG = os.path.join(IMG_DIR, 'ccfs')
    MIL_IMG = os.path.join(IMG_DIR, 'milbourne')


class MovDir:
    """
    Movie directories
    """
    MOV_DIR = os.path.join(BaseDir.BASE_DIR, 'movies')
    NEID_MOV = os.path.join(MOV_DIR, 'neid')
    MERCURY_MOV = os.path.join(MOV_DIR, 'mercury')
    VENUS_MOV = os.path.join(MOV_DIR, 'venus')
    CCF_MOV = os.path.join(MOV_DIR, 'ccfs')


class PlotDir:
    """
    Directories that hold plotting style files
    """
    MPL = os.path.join(BaseDir.BASE_DIR, 'mplstyle')


class Config:
    """
    Config class for usage with CCF Fortran code
    """
    config = {
        'velocity_min': -100,  # km/s minumum velocity at which to compute the CCF
        'velocity_max': 101,  # km/s  minumum velocity at which to compute the CCF
        'velocity_step': 0.25,  # km/s velocity step size at which to compute the CCF
        'qrv': 0,  # km/s systemic velocity of target star
        'mask_name': 'G2_espresso.txt',  # name of mask file containing line list
        'mask_half_width': 0.5,  # km/s width of binary mask (same for all lines)
        'mask_environment': 'air',  # air or vac wavelengths for mask lines
        'ordmin': 57,  # min echelle order number for which to compute CCF (reddest)
        'ordmax': 169,  # max echelle order number for which to compute CCF (bluest)
        'bluest_order': 173,  # bluest order currently extracted (don't change this)
        'resolution': 120000,  # spectral resolution of instrument, lambda / delta_lambda
        'sampling': 5,  # approximate pixel sampling of resolution element (pixels, FWHM)
        'velocity_half_range_to_fit': 20,
        # half width of velocity window for which to fit CCF and derive RV (subset of fill velocity range)
        'clip_edge_pixels': 500  # toss out pixels towards the edge of the array
    }


class Scaling:
    """
    Class that holds scaling coefficients for RV calculations.

    """


class HARPSN(Scaling):
    """
    Class that holds HARPS-N scaling coefficients for RV calculations.

    """

    A = 2.101
    B = 0.9825
    RV0 = 97.08


class NEID(Scaling):
    """
    Class that holds NEID scaling coefficients for RV calculations.

    """

    A = 1.0983
    B = 1.423
    RV0 = -646.076

