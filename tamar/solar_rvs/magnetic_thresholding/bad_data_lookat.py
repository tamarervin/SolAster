"""
Tamar Ervin
Date: July 6, 2021

Looking at the bad magnetic data from NEID
dates.
"""
import datetime
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.animation as animation

import astropy.units as u

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

import tamar.tools.utilities as utils
import tamar.tools.coord_funcs as ctfuncs
import tamar.tools.solar_funcs as sfuncs

csv_file_path = '/example_notebooks_and_data/solar_spectra_fits_files/example_neid_file_timestamps.csv'

# AIA image querying parameters
instrument = a.Instrument.aia
wavelength_list = [a.Wavelength(171 * u.angstrom), a.Wavelength(193 * u.angstrom), a.Wavelength(304 * u.angstrom),
                   a.Wavelength(335 * u.angstrom)]
cadence = a.Sample(45*u.second)  # querying cadence
start_date = '2020-12-30T17:00:00'  # start date of query
end_date = '2020-12-30T17:10:00'

# get data
results = Fido.search(a.Time(start_date, end_date), a.Instrument.hmi, a.Physobs.los_magnetic_field, cadence)

# get files
downloaded_files = (Fido.fetch(results))

# sort files by date
file_download = sorted(downloaded_files)

# make map sequence
mmaps = sunpy.map.Map(file_download)

# transformation for magnetic maps
Bobs = []
for i, mmap in enumerate(mmaps):
    mtheta, mphi, mmu = ctfuncs.coord_trans(mmap, R0=mmap.meta['rsun_obs'] * 760e3 / mmap.meta['rsun_ref'],
                                            outside_map_val=np.nan)
    mmu, mmap = ctfuncs.fix_mu(mmu, mmap)

    B, Br = sfuncs.mag_field(mmu, mmap, B_noise=8)

    Bobs.append(B)

    plt.imshow(B)
    plt.title(mmap.meta['date-obs'])
    plt.savefig('/Users/tervin/images/12_30_17/Bobs_' + str(i))



