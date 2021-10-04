#!/usr/bin/env python
"""
Tamar Ervin
June 16, 2021

automated movie of different HMI bands

HMI Bands (A): 171 (gold), 193 (bronze), 304 (red), 211 (purple), 131 (teal), 335 (blue), 094 (green),
1600 (yellow/green), 1700 (pink)
HMI Products: Magnetogram, Intensitygram, Dopplergram

!!! TODO: make OS independent (deal with paths)

"""


import os

import sys

sys.path.append('/Users/tervin/NEID_Solar_analysis')

import os
import datetime

import matplotlib.pyplot as plt
import matplotlib.animation as animation

import astropy.units as u

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

import tamar.tools.utilities as utils
from tamar.tools.settings import MovDir

# file path to save movie
movie_path = os.path.join(MovDir.MERCURY_MOV, '2019_int.mp4')
fps = 2  # frame rate

# AIA image querying parameters
instrument = a.Instrument.aia
wavelength_list = [a.Wavelength(171 * u.angstrom), a.Wavelength(193 * u.angstrom), a.Wavelength(304 * u.angstrom),
                   a.Wavelength(335 * u.angstrom)]
cadence = a.Sample(1*u.hour)  # querying cadence
start_date = '2019-11-11T00:00:00'  # start date of query
end_date = '2019-11-12T00:00:00'  # end date of query

# # read in csv file to get pertinent dates
# dates_list_long = utils.read_csv(csv_file_path)
#
# # lets just choose a couple dates
# dates_list = dates_list_long


##### ----- build video for ONE wavelength
# if isinstance(dates_list[0], str):
#     datetimes = [datetime.datetime.strptime(date[0], '%Y-%m-%dT%H:%M:%S.%f') for date in dates_list]
# else:
#     datetimes = dates_list
#
# results = []
# for ind, datetime_object in enumerate(datetimes):
#     # if ind % 10 == 0:
#     # pull image within specified time range
#     results.append(Fido.search(a.Time(str(datetime_object - time_range), str(datetime_object + time_range)),
#                                instrument, wavelength_list[1]))
#
# downloaded_files = []
# for ind, datetime_object in enumerate(datetimes):
#     # if ind % 10 == 0:
#     # add file to list
#     downloaded_files.append(Fido.fetch(results[ind]))
#     # file = Fido.fetch(results[ind])


results = Fido.search(a.Time(start_date, end_date), a.Instrument.hmi, a.Physobs.intensity, cadence)
# results = Fido.search(a.Time(start_date, end_date), a.Instrument.hmi, a.Physobs.los_magnetic_field, cadence)
#
# results = Fido.search(a.Time(start_date, end_date), a.Instrument.aia, a.Wavelength(211 * u.angstrom), cadence)


downloaded_files = (Fido.fetch(results))

# check for good results
# remove unusable file types
good_files = []
for file in downloaded_files:
    name, extension = os.path.splitext(file)
    if extension == '.fits':
        good_files.append(file)

# sort files by date
file_download = sorted(good_files)

# make map sequence
file_download = file_download[0:58] + file_download[59:]
aia = sunpy.map.Map(file_download)
aia = [ai.rotate(order=3) for ai in aia]

# limb brightening
import tamar.tools.solar_funcs as sfuncs
import tamar.tools.lbc_funcs as lbfuncs
from sunpy.coordinates import frames

# map_int_cor = []
# for imap in aia:
#     Lij = lbfuncs.limb_polynomial(imap)
#
#     # calculate corrected data
#     Iflat = imap.data / Lij
#
#     # corrected intensity maps
#     map_int_cor.append(sfuncs.corrected_map(Iflat, imap, map_type='Corrected-Intensitygram'))

fig, ax = plt.subplots()
for m in aia:
    m.plot_settings['cmap'] = plt.get_cmap('hinodesotintensity')
    # m.plot_settings['cmap'] = plt.get_cmap('hmimag')
# image half down sequence to get better scaling
plot_obj = aia[len(aia) // 2].plot()


def animate(i):
    ax.set_title("HMI %s %s" % (str(aia[i].meta['content']),
                                aia[i].meta['t_obs']))
    plot_obj.set_data(aia[i].data)
    return (plot_obj,)

# def animate(i):
#     ax.set_title("AIA %s %s" % (str(aia[i].meta['wavelnth']),
#                                 aia[i].meta['t_obs']))
#     plot_obj.set_data(aia[i].data)
#     return (plot_obj,)


anim = animation.FuncAnimation(fig, animate, init_func=None,
                               frames=len(aia), interval=500, blit=True)
Writer = animation.writers['ffmpeg']
writer = Writer(fps, metadata=dict(artist='Me'), bitrate=1800)
anim.save(movie_path, writer=writer)
plt.close(fig)

# for i, img in enumerate(aia):
#     img.plot_settings['norm'] = plt.Normalize(-100, 100)
#     img.plot()
#     plt.savefig('/Users/tervin/images/12_30_17/img_' + str(i) + '.png')

# ffmpeg -r 1 -i img_%d.png -c:v libx264 -vf "fps=25,format=yuv420p" 12_30_17_mag.mp4

##### ----- build video for multiple wavelengths

# # download images for all wavelengths
# results = []
# for wave_ind, wavelength in enumerate(wavelength_list[1:]):
#     results_wave = []
#     for date_ind, datetime_object in enumerate(datetimes):
#         # pull image within specified time range
#         results_wave.append(Fido.search(a.Time(str(datetime_object - time_range), str(datetime_object + time_range)),
#                                         instrument, wavelength))
#     results.append(results_wave)
#
# downloaded_files = []
# for ind, datetime_object in enumerate(datetimes):
#         # add file to list
#         downloaded_files.append(Fido.fetch(results[ind]))
#
#
# # make map sequences for each wavelength
# img = []
# for im in maps_list:
#     fig, axs = plt.subplots(2, 2)
#     axs[0, 0].plot(maps_list[2])
#     axs[0, 1].plot(maps_list[1])
#     axs[1, 0].plot(maps_list[2])
#     axs[1, 1].plot(maps_list[3])
#     # img.append([im.plot()])
