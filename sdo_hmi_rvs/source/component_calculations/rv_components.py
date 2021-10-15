"""
Tamar Ervin
Date: June 22, 2021

calculation of solar RV components

"""

import os
import sys
import numpy as np

import matplotlib.pyplot as plt

import astropy.units as u

import sunpy.map
from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.coordinates import frames

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../'))
import sdo_hmi_rvs.tools.calculation_funcs as sfuncs
import sdo_hmi_rvs.tools.lbc_funcs as lbfuncs
import sdo_hmi_rvs.tools.coord_funcs as ctfuncs

# get hmi data products
cadence = a.Sample(24*u.hour)  # querying cadence
start_date = '2021-01-01T12:00:00'  # start date of query
end_date = '2021-01-02T12:00:00'
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
        hmi_vel.append(map_obj.rotate(order=3))
    elif map_obj.meta['content'] == 'MAGNETOGRAM':
        hmi_mag.append(map_obj.rotate(order=3))
    elif map_obj.meta['content'] == 'CONTINUUM INTENSITY':
        hmi_int.append(map_obj.rotate(order=3))

# transformation for velocity maps
vcoords = [ctfuncs.coord_trans(smap, R0=smap.meta['rsun_obs'] * 760e3 / smap.meta['rsun_ref'], outside_map_val=np.nan) for smap in hmi_vel]

# transformation for intensity maps
icoords = [ctfuncs.coord_trans(smap, R0=smap.meta['rsun_obs'] * 760e3 / smap.meta['rsun_ref'], outside_map_val=np.nan) for smap in hmi_int]

# transformation for magnetic maps
mcoords = [ctfuncs.coord_trans(smap, R0=smap.meta['rsun_obs'] * 760e3 / smap.meta['rsun_ref'], outside_map_val=np.nan) for smap in hmi_mag]

# calculate relative positions
delta = [sfuncs.rel_positions(vcoords[i][0], vcoords[i][1], hmi_vel[i]) for i in range(0, len(hmi_vel))]

# calculate spacecraft velocity
vsc = [sfuncs.spacecraft_vel(delta[i][0], delta[i][1], delta[i][2], delta[i][3], hmi_vel[i]) for i in range(0, len(hmi_vel))]

# optimized solar rotation parameters
a1 = 14.713
a2 = -2.396
a3 = -1.787
a_parameters = [a1, a2, a3]

# calculation of solar rotation velocity
vrot = [sfuncs.solar_rot_vel(a_parameters, vcoords[i][0], delta[i][0], delta[i][1], delta[i][2], delta[i][3], hmi_vel[i]) for i in range(0, len(hmi_vel))]

# calculate corrected velocity
corrected_vel = [hmi_vel[i].data - vrot[i] - vsc[i] for i in range(0, len(hmi_vel))]

# corrected velocity maps
map_vel_cor = [sfuncs.corrected_map(corrected_vel[i], hmi_vel[i], map_type='Corrected-Dopplergram', frame=frames.HeliographicCarrington) for i in range(0, len(hmi_vel))]

# limb brightening
Lij = [lbfuncs.limb_polynomial(hmi_int[i]) for i in range(0, len(hmi_int))]

# calculate corrected data
Iflat = [hmi_int[i].data/Lij[i] for i in range(0, len(hmi_int))]

# corrected intensity maps
map_int_cor = [sfuncs.corrected_map(Iflat[i], hmi_int[i], map_type='Corrected-Intensitygram', frame=frames.HeliographicCarrington) for i in range(0, len(hmi_int))]

# magnetic noise level
B_noise = 8

# calculate unsigned field strength
corrected_mag = [sfuncs.mag_field(mcoords[i][2], hmi_mag[i], B_noise) for i in range(0, len(hmi_mag))]

# corrected observed magnetic data map
map_mag_obs = [sfuncs.corrected_map(corrected_mag[i][0], hmi_mag[i], map_type='Corrected-Magnetogram', frame=frames.HeliographicCarrington) for i in range(0, len(hmi_mag))]

# radial magnetic data map
map_mag_cor = [sfuncs.corrected_map(corrected_mag[i][1], hmi_mag[i], map_type='Corrected-Magnetogram', frame=frames.HeliographicCarrington) for i in range(0, len(hmi_mag))]

### calculate magnetic threshold
# magnetic threshold value (G) from Yeo et. al. 2013
Br_cutoff = 24
# mu cutoff value
mu_cutoff = 0.3
# calculate magnetic threshold
mag_thresh = [sfuncs.mag_thresh(mcoords[i][2], map_mag_cor[i], Br_cutoff=Br_cutoff, mu_cutoff=mu_cutoff) for i in range(0, len(hmi_mag))]

# thresholded magnetic maps
map_mag_thresh = [sfuncs.corrected_map(mag_thresh[i][0], hmi_mag[i], map_type='Magnetic-Threshold', frame=frames.HeliographicCarrington) for i in range(0, len(hmi_mag))]

# calculate intensity threshold
int_thresh = [sfuncs.int_thresh(map_int_cor[i], map_mag_cor[i], mag_thresh[i][1]) for i in range(0, len(hmi_int))]

# create threshold array
thresh_arr = [sfuncs.thresh_map(int_thresh[i][2], int_thresh[i][3]) for i in range(0, len(hmi_mag))]

# full threshold maps
map_full_thresh = [sfuncs.corrected_map(thresh_arr[i], hmi_mag[i], map_type='Threshold', frame=frames.HeliographicCarrington) for i in range(0, len(hmi_int))]

### velocity contribution due to convective motion of quiet-Sun
v_quiet = [sfuncs.v_quiet(map_vel_cor[i], hmi_int[i], mag_thresh[i][1]) for i in range(0, len(hmi_mag))]

### velocity contribution due to rotational Doppler imbalance of active regions (faculae/sunspots)
# create array of active weights (flipped from mag weights)
active_weights = [sfuncs.active_weights(mag_thresh[i][1]) for i in range(0, len(hmi_int))]

# calculate photospheric velocity
v_phot = [sfuncs.v_phot(mag_thresh[i][1], active_weights[i], Lij[i], vsc[i], hmi_int[i], hmi_vel[i]) for i in range(0, len(hmi_int))]

### velocity contribution due to suppression of convective blueshift by active regions
# calculate disc-averaged velocity
v_disc = [sfuncs.v_disc(map_vel_cor[i], hmi_int[i]) for i in range(0, len(hmi_int))]

# calculate convective velocity
v_conv = [v_disc[i] - v_quiet[i] - v_phot[i] for i in range(0, len(hmi_mag))]

# calculate filling factor
filling_factor = [sfuncs.filling_factor(hmi_mag[i], active_weights[i]) for i in range(0, len(hmi_mag))]

### unsigned magnetic flux
# unsigned observed flux
unsigned_obs_flux = [sfuncs.unsigned_flux(map_mag_obs[i], hmi_int[i]) for i in range(0, len(hmi_mag))]

# unsigned radial flux
unsigned_rad_flux = [sfuncs.unsigned_flux(map_mag_obs[i], hmi_int[i]) for i in range(0, len(hmi_mag))]

#### ----- plot some stuff ----- ####

#### ----- magnetic maps
mag_map = hmi_mag[0]
mag_map.plot_settings['norm'] = plt.Normalize(-100, 100)
mag_map.plot()
plt.colorbar()
plt.show()

mag_map = map_mag_cor[0]
mag_map.plot_settings['norm'] = plt.Normalize(-100, 100)
mag_map.plot()
plt.colorbar()
plt.show()

mag_map.peek()

#### ----- velocity maps
vel_map = hmi_vel[0]
vel_map.plot_settings['norm'] = plt.Normalize(-5500, 5500)
vel_map.plot()
plt.colorbar()
plt.show()

vel_map1 = map_vel_cor[0]
vel_map1.plot()
plt.colorbar()
plt.show()

#### ----- intensity maps
int_map = hmi_int[0]
cont_map = plt.get_cmap('hinodesotintensity')
int_map.plot_settings['cmap'] = cont_map
int_map.plot()
plt.colorbar()
plt.show()

int_map = map_int_cor[0]
cont_map = plt.get_cmap('hinodesotintensity')
int_map.plot_settings['cmap'] = cont_map
int_map.plot()
plt.colorbar()
plt.show()

hmi_int[0].peek()
map_int_cor[0].peek()
#### ----- faculae to sunspot comparison
plt.imshow(map_full_thresh[0].data, cmap=plt.get_cmap('bwr'))
plt.colorbar()
plt.title("Map showing faculae (blue) and sunspots (red)")

#### ----- compare velocity with flattened intensity

# index of image you want to use
i = 0

# setup plot
plt.figure("Velocity v. Intensity")
plt.title("LOS Velocity as function of Flattened Intensity: \n" +
          str(hmi_mag[i].meta['date-obs']))
plt.xlabel("Line-of-Sight Velocity")
plt.ylabel("Flattened Continuum Intensity")

# plot non-active pixels in black
x = hmi_vel[i].data
y = Iflat[i]
plt.scatter(x*mag_thresh[i][1], y*mag_thresh[i][1], color='k', label='Non-Active')

# plot active pixels in green
plt.scatter(x*active_weights[i], y*active_weights[i], color='green', label='Active')

# horizontal intensity threshold line
plt.axhline(int_thresh[i][4], xmin=x.min(), xmax=x.max(), color='purple', linestyle='dashed')

# legend
plt.legend()
plt.show()

#### ----- compare observed magnetic field to flattened intensity

# index of image you want to use
i = 0

# get x and y data to plot
x = np.abs(map_mag_obs[i].data)
y = Iflat[i]

# setup plot
plt.figure("Magnetic Field v. Intensity")
plt.title("Observed Magnetic Field as function of Flattened Intensity: \n" +
          str(hmi_mag[i].meta['date-obs']))
plt.xlabel("Observed Magnetic Field")
plt.ylabel("Flattened Continuum Intensity")

# plot non-active pixels in black
plt.scatter(x*mag_thresh[i][1], y*mag_thresh[i][1], color='k', label='Non-Active')

# plot active pixels in green
plt.scatter(x*active_weights[i], y*active_weights[i], color='green', label='Active')

# legend
plt.legend()
plt.show()


