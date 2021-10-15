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
    BASE_DIR = ''


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

