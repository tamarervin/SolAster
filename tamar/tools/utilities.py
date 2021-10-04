"""
Tamar Ervin
Date: June 24, 2021

Utility functions for reading and writing data

"""

import csv
import datetime

import numpy as np
import pandas as pd

from astropy.time import Time


def read_csv(csv_file_path):
    """
    function to read csv file and return list of dates

    Parameters
    ----------
    csv_file_path: path to csv file (string)

    Returns
    -------
    csv_list: list (strings) of elements in csv file

    """

    with open(csv_file_path, newline='') as f:
        reader = csv.reader(f)
        csv_list = list(reader)

    return csv_list


def append_list_as_row(file_name, list_of_elem):
    """
    function to add row to csv file

    Parameters
    ----------
    file_name: path to csv file
    list_of_elem: elements as a list to add to file

    Returns
    -------

    """
    # Open file in append mode
    with open(file_name, 'a+', newline='') as write_obj:
        # Create a writer object from csv module
        csv_writer = csv.writer(write_obj)
        # Add contents of list as last row in the csv file
        csv_writer.writerow(list_of_elem)

    return None


def get_dates(date):
    """
    function to convert dates from either JD, string, or datetime
    to a Sunpy usable date form

    Parameters
    ----------
    date: date in any form (JD, string, datetime)

    Returns
    -------
    date_str: UT datetime as string
    date_obj: UT datetime object

    """
    if isinstance(date, str):
        date_str = date
        date_obj = datetime.datetime.strptime(date_str, '%Y-%m-%dT%H:%M:%S')
        date_jd = Time(date_str)
        date_jd.format = 'jd'
    elif isinstance(date, float):
        t = Time(date, format='jd')
        date_obj = t.datetime
        date_str = date_obj.strftime('%Y-%m-%dT%H:%M:%S')
        date_jd = date
    else:
        date_obj = date
        date_str = date_obj.strftime('%Y-%m-%dT%H:%M:%S')
        date_jd = Time(date_str)
        date_jd.format = 'jd'

    return date_str, date_obj, date_jd


def get_neid_struct_dates(date):
    """
    function to convert dates to the format used by the
    NEID fits files directory structure

    Parameters
    ----------
    date: date in any form (JD, string, datetime)

    Returns
    -------

    """

    if isinstance(date, str):
        # date_str = date[0:4] + date[5:7] + date[8:10]
        date_str = date
        date_obj = datetime.datetime.strptime(date_str, '%Y%m%d')
        date_str = date_obj.strftime('%Y-%m-%dT%H:%M:%S')
        date_jd = Time(date_str)
        date_jd.format = 'jd'
    elif isinstance(date, float):
        t = Time(date, format='jd')
        date_obj = t.datetime
        date_str = date_obj.strftime('%Y%m%d')
        date_jd = date
    else:
        date_obj = date
        date_str = date_obj.strftime('%Y%m%d')
        date_jd = Time(date_str)
        date_jd.format = 'jd'

    return date_str, date_obj, date_jd


def read_sdo_csv(csv_file):
    """
    function to read in csv file with sdo calculations
    and return metrics

    Parameters
    ----------
    csv_file: path to csv file to read

    Returns
    -------
    list of SDO derived metrics

    """

    # create pandas dataframe
    component_df = pd.read_csv(csv_file)

    # get dates list
    date_jd = component_df.date_jd.values
    inds = np.argsort(date_jd)
    rv_sun = component_df.rv_sun.values[inds]
    rv_error = component_df.rv_error.values[inds]

    non_nan = np.logical_not(np.isnan(rv_sun))

    rv_med = np.nanmedian(np.abs(rv_sun))

    good_sun = np.logical_and(np.abs(rv_sun) > rv_med - 2, np.abs(rv_sun) < rv_med + 2)
    good_error = np.logical_and(np.abs(rv_error) < .4, np.abs(rv_error) < 0.37)
    good = np.logical_and(good_sun, good_error)

    good_rvs = np.logical_and(good, non_nan)

    # get velocities lists
    v_phot = component_df.v_phot.values[inds][good_rvs]
    v_conv = component_df.v_conv.values[inds][good_rvs]
    rv_model = component_df.rv_model.values[inds][good_rvs]
    rv_sun = rv_sun[good_rvs]
    rv_error = rv_error[good_rvs]

    # get magnetic observables
    f = component_df.f.values[inds][good_rvs]
    Bobs = component_df.Bobs.values[inds][good_rvs]
    f_bright = component_df.f_bright.values[inds][good_rvs]
    f_spot = component_df.f_spot.values[inds][good_rvs]

    # dates
    date_jd = date_jd[inds][good_rvs]
    dates = component_df.date_obs.values[inds][good_rvs]

    return date_jd, dates, v_phot, v_conv, rv_model, rv_sun, rv_error, f, Bobs, f_bright, f_spot


def read_ccf_csv(csv_file):
    """
    function to read in csv file with ccf calculations
    and return metrics

    Parameters
    ----------
    csv_file: path to csv file to read

    Returns
    -------
    list of CCF derived metrics

    """
    fpickle = pd.read_pickle(csv_file)
    dates = fpickle.dates.values
    jd_dates = fpickle.jd_dates.values

    inds = np.argsort(jd_dates)
    rv_sun = fpickle.rv_sun.values[inds] * 1e3
    rv_error = fpickle.rv_error.values[inds] * 1e3

    non_nan = np.logical_not(np.isnan(rv_sun))

    rv_med = np.nanmedian(np.abs(rv_sun))

    good_sun = np.logical_and(np.abs(rv_sun) > rv_med - 2, np.abs(rv_sun) < rv_med + 2)
    good_error = np.logical_and(np.abs(rv_error) < .4, np.abs(rv_error) < 0.37)
    good = np.logical_and(good_sun, good_error)

    good_rvs = np.logical_and(good, non_nan)

    dates = dates[inds][good_rvs]
    jd_dates = jd_dates[inds][good_rvs]
    rv_sun = rv_sun[good_rvs]
    rv_error = rv_error[good_rvs]
    ccf_list = fpickle.ccf.values[inds][good_rvs]
    rv_model = fpickle.rv_model.values[inds][good_rvs]
    gaussian = fpickle.gaussian.values[inds][good_rvs]
    rv_gauss = fpickle.rv_gauss.values[inds][good_rvs]
    skew = fpickle.median_skew.values[inds][good_rvs]
    int_area = fpickle.integrated_area.values[inds][good_rvs]
    Bobs = fpickle.Bobs.values[inds][good_rvs]
    v_conv = fpickle.v_conv.values[inds][good_rvs]

    # find bad CCFs -- need to figure out why this is...
    bad_ccfs = np.array([np.isnan(ccf[0]) for ccf in ccf_list])
    dates = dates[~bad_ccfs]
    jd_dates = jd_dates[~bad_ccfs]
    ccf_list = ccf_list[~bad_ccfs]
    rv_sun = rv_sun[~bad_ccfs]
    rv_model = rv_model[~bad_ccfs] * 1e3
    rv_error = rv_error[~bad_ccfs]
    gaussian = gaussian[~bad_ccfs]
    rv_gauss = rv_gauss[~bad_ccfs] * 1e3
    skew = skew[~bad_ccfs] * 1e3
    int_area = int_area[~bad_ccfs]
    Bobs = Bobs[~bad_ccfs]
    v_conv = v_conv[~bad_ccfs]

    return dates, jd_dates, ccf_list, rv_sun, rv_model, rv_error, gaussian, rv_gauss, skew, int_area, Bobs, v_conv

