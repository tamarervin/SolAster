"""
Tamar Ervin
Date: June 17, 2021

outline for step 1 in calculation of solar RVs
based on method from Haywood et. al. 2016
- coordinate transformation

"""

import numpy as np

import astropy.units as u

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

import tamar.tools.solar_funcs as sfuncs

# get magnetograms
cadence = a.Sample(24*u.hour)  # querying cadence
start_date = '2021-01-01T12:00:00'  # start date of query
end_date = '2021-01-05T12:00:00'

result = Fido.search(a.Time(start_date, end_date),
                     a.Instrument.hmi, a.Physobs.los_magnetic_field, cadence)


# download results
file_download = Fido.fetch(result)
# convert to map sequence
map_seq = sunpy.map.Map(sorted(file_download))

# convert to map sequence
map_seq = sunpy.map.Map(sorted(file_download))
hmi_vel, hmi_mag, hmi_int \
    = [], [], []
for i, map_obj in enumerate(map_seq):
    if map_obj.meta['content'] == 'DOPPLERGRAM':
        hmi_vel.append(map_obj)
    elif map_obj.meta['content'] == 'MAGNETOGRAM':
        hmi_mag.append(map_obj)
    elif map_obj.meta['content'] == 'CONTINUUM INTENSITY':
        hmi_int.append(map_obj)

# choose one map for starters
vel_map = hmi_vel[0].rotate(order=3)

# rsun = np.float32(vel_map.rsun_meters)
rsun = vel_map.meta['rsun_obs'] * 760e3 / vel_map.meta['rsun_ref']
x, y = sfuncs.get_scales_from_map(vel_map)
theta, phi, mu = sfuncs.get_coordinates(x, y, vel_map, R0=rsun, outside_map_val=np.nan)

# # downsample maps to take up less memory
# maps = [m.resample((1024, 1024)*u.pix) for m in map_seq]
#
# # convert HMI magnetogram from pixel coordinates to heliographic coordinates
# hmi_map = maps[0]
#
# # make new header file for updated coordinate frames
# shape_out = (1024, 1024)
# header = sunpy.map.make_fitswcs_header(shape_out,
#                                        SkyCoord(0, 0, unit=u.deg,
#                                                 frame=frames.HeliographicCarrington,
#                                                 obstime=maps[0].date, observer=maps[0].observer_coordinate),
#                                        scale=[180 / shape_out[0],
#                                               360 / shape_out[1]] * u.deg / u.pix,
#                                        wavelength=int(maps[0].meta['wavelnth']) * u.AA,
#                                        projection_code="CAR")
# out_wcs = WCS(header)
# # coordinate transformation
# coordinates = tuple(map(sunpy.map.all_coordinates_from_map, maps))
# weights = [coord.transform_to(frames.HeliographicCarrington) for coord in coordinates]
# # smooth out the transformation
# weights = [(w / np.nanmax(w)) ** 3 for w in weights]
# for w in weights:
#     w[np.isnan(w)] = 0
#
#
# hmi_coord_trans = sunpy.coordinates.transformations.hgs_to_hgc(hmi_map.observer_coordinate, frames.HeliographicCarrington(observer=hmi_map.observer_coordinate))
#
# ##### ugh do this by coordinate?? so tedious
# hpc1 = SkyCoord(0*u.arcsec, 0*u.arcsec, frame=hmi_map.coordinate_frame)
# hpc2 = [hpc1.transform_to(frames.HeliographicCarrington) for hpc1 in hmi_map]
#
# # reproject
# array, _ = reproject_and_coadd(maps, out_wcs, shape_out,
#                                input_weights=weights,
#                                reproject_function=reproject_interp,
#                                match_background=True,
#                                background_reference=0)
# array, footprint = reproject_and_coadd(maps, out_wcs, shape_out,
#                                        reproject_function=reproject_interp)
# reproject_maps = sunpy.map.Map((array, header))

### old tests byeeee
# map_seq_trans = [map_obj.coordinate_frame.transform_to(frames.HeliographicStonyhurst) for map_obj in map_seq]
# coor = map_seq[2].coordinate_frame
# coor.is_transformable_to(frames.HeliographicCarrington)
# coor1 = map_seq[2].coordinate_frame.transform_to(frames.HeliographicCarrington)
#
#
# map_seq[1].coordinate_frame.transform_to(frames.HeliographicStonyhurst)
# sc = map_seq[1].data.pixel_to_data(100 * u.pixel, 100 * u.pixel)
# sc.transform_to(frames.HeliographicCarrington)
#
# import sunpy.map.maputils as maputils
# map_pixels = maputils.all_pixel_indices_from_map(map_seq[1])
# map_coords = maputils.all_coordinates_from_map(map_seq[1])
# map_coords.transform_to(frames.HeliographicCarrington)


####

# sc = SkyCoord(x*u.deg, y*u.deg, frame=frames.HeliographicCarrington, observer=vel_map.observer_coordinate)
#
# # convert coordinates - shaky
# coord = SkyCoord(0*u.arcsec, 0*u.arcsec, obstime=vel_map.date, observer=vel_map.observer_coordinate, frame=frames.HeliographicCarrington)
# header = sunpy.map.make_fitswcs_header(vel_map.data, coord)
# hmi_hgc_map = sunpy.map.Map(vel_map.data, header)
# hmi_hgc_map.peek()
# x, y = get_scales_from_map(hmi_hgc_map)
#
# test = hmi_map.data_to_pixel(hmi_map.data, 1)
# # get coordinate grid
# botw = hmi_hgc_map.bottom_left_coord.lat
# botn = hmi_hgc_map.bottom_left_coord.lon
# botr = hmi_hgc_map.bottom_left_coord.radius
#
# bot = hmi_hgc_map.bottom_left_coord.T
#
# # round 2 conversion
# bot_left = hmi_map.bottom_left_coord
# bot_hgc = bot_left.transform_to(frame=frames.HeliographicCarrington)
# top_right = hmi_map.top_right_coord
# top_hgc = top_right.transform_to(frame=frames.HeliographicCarrington)
#
# hgc_lat = np.linspace(hmi_hgc_map.bottom_left_coord.lat, hmi_hgc_map.top_right_coord.lat, 4096)
# hgc_lon = np.linspace(hmi_hgc_map.bottom_left_coord.lon, hmi_hgc_map.top_right_coord.lon, 4096)
# hgc_rad = np.linspace(hmi_hgc_map.bottom_left_coord.rad, hmi_hgc_map.top_right_coord.rad, 4096)
#
# sc = SkyCoord(hgc_lat, hgc_lon, hgc_rad, frame=frames.HeliographicCarrington)
#
# botx = bot_left.x
# boty = bot_left.y
# botz = bot_left.z
# topx = top_right.x
# topy = top_right.y
# topz = top_right.z
#
# x = np.linspace(botx, topx, 4096) * u.deg
# y = np.linspace(boty, topy, 4096) * u.deg
# z = np.linspace(botz, topz, 4096) * u.deg
#
# x1 = np.linspace(hmi_map.bottom_left_coord.Tx, hmi_map.top_right_coord.Tx, 4096)
# y1 = np.linspace(hmi_map.bottom_left_coord.Ty, hmi_map.top_right_coord.Ty, 4096)
#
# sc_list = [SkyCoord(xi, yi, frame=frames.Helioprojective, observer=hmi_map.observer_coordinate) for xi, yi in zip(x1, y1)]
# sc_trans = [sc.transform_to(frame=frames.HeliographicCarrington) for sc in sc_list]
#
# sc = SkyCoord(x1, y1, frame=frames.Helioprojective, observer=hmi_map.observer_coordinate)
# sc1 = sc.transform_to(frame=frames.HeliographicCarrington)

# load fits info
# from astropy.io import fits
# def load_wcs_from_file(map):
#
#     # Parse the WCS keywords in the primary HDU
#     w = WCS(map.fits_header)
#
#     # Print out all of the settings that were parsed from the header
#     # w.wcs.print_contents()
#
#     # Three pixel coordinates of interest.
#     # Note we've silently assumed an NAXIS=2 image here.
#     # The pixel coordinates are pairs of [X, Y].
#     # The "origin" argument indicates whether the input coordinates
#     # are 0-based (as in Numpy arrays) or
#     # 1-based (as in the FITS convention, for example coordinates
#     # coming from DS9).
#     idx = np.indices((4096, 4096))
#     idx_row = idx[0].flatten()
#     idx_col = idx[1].flatten()
#     pixcrd = np.array(idx_row, idx_col)
#
#     # Convert pixel coordinates to world coordinates
#     # The second argument is "origin" -- in this case we're declaring we
#     # have 0-based (Numpy-like) coordinates.
#     world = w.wcs_pix2world(pixcrd, 0)
#
#     # Convert the same coordinates back to pixel coordinates.
#     pixcrd2 = w.wcs_world2pix(world, 0)
#
#     # These should be the same as the original pixel coordinates, modulo
#     # some floating-point error.
#     assert np.max(np.abs(pixcrd - pixcrd2)) < 1e-6
#
#     # The example below illustrates the use of "origin" to convert between
#     # 0- and 1- based coordinates when executing the forward and backward
#     # WCS transform.
#     x = 0
#     y = 0
#     origin = 0
#     assert (w.wcs_pix2world(x, y, origin) ==
#             w.wcs_pix2world(x + 1, y + 1, origin + 1))
#
#     return world
#
# test = load_wcs_from_file(hmi_map)
