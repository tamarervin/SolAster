"""
Tamar Ervin
June 17, 2021

automated movie of different HMI bands

HMI Bands (A): 171 (gold), 193 (bronze), 304 (red), 211 (purple), 131 (teal), 335 (blue), 094 (green),
1600 (yellow/green), 1700 (pink)
HMI Products: Magnetogram, Intensitygram, Dopplergram

!!! TODO: make OS independent (deal with paths)

"""

#### THIS TIMESTAMP WORKS TO MAKE MOVIE AND SUBPLOTS (June 17, 2021 @ 11:58)

import matplotlib.pyplot as plt
import matplotlib.animation as animation

import astropy.units as u

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

# file path to save movie
movie_path = '/Users/tervin/NEID_Solar_analysis/tamar/aia_waves.mp4'
fps = 2  # frame rate

# AIA image querying parameters
instrument = a.Instrument.aia
wavelength_list = [a.Wavelength(171 * u.angstrom), a.Wavelength(193 * u.angstrom), a.Wavelength(304 * u.angstrom),
                   a.Wavelength(335 * u.angstrom)]
query_range = a.Sample(24*u.hour)
start_date = '2021-01-01T12:00:00'
end_date = '2021-01-10T12:00:00'

# get results for each wavelength
results = []
downloaded_files = []
for wavelength in wavelength_list:
    result = Fido.search(a.Time(start_date, end_date), instrument, wavelength, query_range)
    file_download = Fido.fetch(result)
    results.append(results)
    downloaded_files.append(file_download)


# make movie
aia = []
for files in downloaded_files:
    aia_file = sunpy.map.Map(files[0:11])
    aia.append(aia_file)

fig, axs = plt.subplots(2, 2)


def animate(l):
    fig.suptitle("AIA %s" % (str(aia[0][l].meta['t_obs'])))
    axs[0, 0].set_title("AIA %s" % (str(aia[0][l].meta['wavelnth'])))
    axs[0, 1].set_title("AIA %s" % (str(aia[1][l].meta['wavelnth'])))
    axs[1, 0].set_title("AIA %s" % (str(aia[2][l].meta['wavelnth'])))
    axs[1, 1].set_title("AIA %s" % (str(aia[3][l].meta['wavelnth'])))
    # for j in range(0, 2):
    #     for k in range(0, 2):
    #         axs[j, k].set_xticks([])
    #         axs[j, k].set_yticks([])
    plot_obj0 = aia[0][l].plot(axes=axs[0, 0])
    plot_obj1 = aia[1][l].plot(axes=axs[0, 1])
    plot_obj2 = aia[2][l].plot(axes=axs[1, 0])
    plot_obj3 = aia[3][l].plot(axes=axs[1, 1])
    # plot_obj0.set_data(aia[0][l].data)
    # plot_obj1.set_data(aia[1][l].data)
    # plot_obj2.set_data(aia[2][l].data)
    # plot_obj3.set_data(aia[3][l].data)

    return plot_obj0, plot_obj1, plot_obj2, plot_obj3


anim = animation.FuncAnimation(fig, animate, init_func=None, frames=len(aia[0]), interval=100, blit=True)
Writer = animation.writers['ffmpeg']
writer = Writer(fps, metadata=dict(artist='Me'), bitrate=1800)

anim.save(movie_path, writer=writer)
plt.close(fig)


### making subplots of sunpy maps
i = 0
fig, axs = plt.subplots(2, 2)
fig.suptitle("AIA %s" % (str(aia[0][i].meta['t_obs'])))
axs[0, 0].set_title("AIA %s" % (str(aia[0][i].meta['wavelnth'])))
axs[0, 1].set_title("AIA %s" % (str(aia[1][i].meta['wavelnth'])))
axs[1, 0].set_title("AIA %s" % (str(aia[2][i].meta['wavelnth'])))
axs[1, 1].set_title("AIA %s" % (str(aia[3][i].meta['wavelnth'])))
for j in range(0, 2):
    for k in range(0, 2):
        axs[j, k].set_xticks([])
        axs[j, k].set_yticks([])
aia[0][i].plot(axes=axs[0, 0], annotate=False)
aia[1][i].plot(axes=axs[0, 1], annotate=False)
aia[2][i].plot(axes=axs[1, 0], annotate=False)
aia[3][i].plot(axes=axs[1, 1], annotate=False)

plt.show()


# rando testing
# def animate(i):
#     fig.suptitle("AIA %s" % (str(aia[0][i].meta['t_obs'])))
#     axs0.set_title("AIA %s" % (str(aia[0][i].meta['wavelnth'])))
#     axs1.set_title("AIA %s" % (str(aia[1][i].meta['wavelnth'])))
#     axs2.set_title("AIA %s" % (str(aia[2][i].meta['wavelnth'])))
#     axs3.set_title("AIA %s" % (str(aia[3][i].meta['wavelnth'])))
#     # plot_obj0 = aia[0][i].plot(axes=axs0)
#     # plot_obj0.set_data(aia[0][i].data)
#     axs0.imshow(aia[0][i].data)
#     axs1.imshow(aia[1][i].data)
#     axs2.imshow(aia[2][i].data)
#     axs3.imshow(aia[3][i].data)
#     # plot_obj1 = aia[1][i].plot(axes=axs1)
#     # plot_obj1.set_data(aia[1][i].data)
#     # plot_obj2 = aia[2][i].plot(axes=axs2)
#     # plot_obj2.set_data(aia[2][i].data)
#     # plot_obj3 = aia[3][i].plot(axes=axs3)
#     # plot_obj3.set_data(aia[3][i].data)
#     return axs0, axs1, axs2, axs3


# import astropy
# from sunpy.visualization.animator import ArrayAnimatorWCS
# map_sequence = sunpy.map.Map(downloaded_files[0], sequence=True)
# sequence_array = map_sequence.as_array()
# wcs = astropy.wcs.WCS(naxis=3)
# wcs_anim = ArrayAnimatorWCS(sequence_array, wcs, [0, 'x', 'y'], vmax=1000)
# plt.show()

