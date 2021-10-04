"""
Tamar Ervin
Date: June 21, 2021

flatten continuum intensity images to correct for limb-darkening
use fifth order polynomial and constants from Allen 1973

IDL code here: https://hesperia.gsfc.nasa.gov/ssw/gen/idl/solar/darklimb_correct.pro
seems there is currently no python version...use stuff from PSI??!
"""

import astropy.units as u

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a

import tamar.tools.coord_funcs as ctfuncs


# get magnetograms
cadence = a.Sample(24*u.hour)  # querying cadence
start_date = '2021-01-01T12:00:00'  # start date of query
end_date = '2021-01-05T12:00:00'
physobs_list = [a.Physobs.los_velocity, a.Physobs.los_magnetic_field, a.Physobs.intensity]

result = Fido.search(a.Time(start_date, end_date),
                     a.Instrument.hmi, physobs_list[0] | physobs_list[1] | physobs_list[2], cadence)


# download results
file_download = Fido.fetch(result)
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

int_map = hmi_int[0].rotate(order=3)

# calculate mu value (from org_vel_corrections.py)
x, y = ctfuncs.get_scales_from_map(int_map)
theta, phi, mu = ctfuncs.get_coordinates(x, y, int_map, R0=1.01, outside_map_val=0)


# TODO: calculate limb brightening correction by whatever method
Iij = int_map.data
Lij = ""  # correction polynomial
Iflat = Iij/Lij