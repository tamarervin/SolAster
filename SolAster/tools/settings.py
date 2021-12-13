"""
Tamar Ervin
Date: June 14, 2021

Settings file for paths
Update to work with your system.
"""

import os


class BaseDir:
    """
    Base directory.
    Update this if you want to save the files outside the repo and follow same file structure and all will be well :)
    """
    BASE_DIR = os.path.realpath('../products/')


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


class PlotDir:
    """
    Directories that hold plotting style files
    """
    # MPL = os.path.join(BaseDir.BASE_DIR, 'mplstyle')
    MPL = '../products/mplystyle/'


class Parameters:
    """
    Class to hold parameters for RV pipeline calculation.
    May be updates by the user.
    """
    # minimum mu value for 'good' pixels
    mu_cutoff = 0.3

    # coefficients for differential rotation from Snodgrass & Ulrich (1990).
    a1 = 14.713
    a2 = -2.396
    a3 = -1.787

    # magnetic noise level used to correct for unsigned magnetic field strength
    B_noise = 8

    # magnetic threshold value (G) from Yeo et. al. 2013 used for identifying magnetic regions
    Br_cutoff = 24

    # area threshold value between large and small regions (in uHem) based on Milbourne et al. (2019)
    athresh = 20


class Scaling:
    """
    Class that holds scaling coefficients for RV calculations.
    """
    pass


class HARPSN(Scaling):
    """
    Class that holds HARPS-N scaling coefficients for RV calculations.
    Fit using linear regression on HARPS-N data from 2015.
    """

    A = 2.101
    B = 0.9825
    RV0 = 97.08


class NEID(Scaling):
    """
    Class that holds NEID scaling coefficients for RV calculations.
     Fit using linear regression on NEID data from Dec 2020 - June 2021.
    """

    A = 1.0983
    B = 1.423
    RV0 = -646.076

