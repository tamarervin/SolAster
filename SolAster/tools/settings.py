"""
Tamar Ervin
Date: June 14, 2021

Settings file for paths
Update to work with your system.
"""

import os
import datetime


class Inputs:
    """
    Class to hold user specified inputs to run examples.
    See README or documentation site for additional information.

    """

    # name of csv file to store calculations
    csv_name = 'paper_test'

    # name of instrument to use for calculation of RV model
    # choose either 'NEID' or 'HARPS-N'
    inst = 'NEID'

    # querying cadence in seconds
    cadence = 24 * 60 * 60

    # start date for calculations
    start_date = datetime.datetime(2021, 2, 7, 0, 00, 0)

    # end date for calculations
    end_date = datetime.datetime(2021, 2, 10, 0, 00, 0)

    # True if outputting diagnostic plots
    diagnostic_plots = True
    # path to save diagnostic figure or none
    save_fig = None


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
    CALC = os.path.join(CSV_DIR, 'calcs')

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

