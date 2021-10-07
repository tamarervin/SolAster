# Data Corrections

For the calculation of independent RV's, we used the three main
HMI data products: Dopplergrams, Intensitygrams, Magnetograms. These data products
must be corrected for inherent issues such as relative velocityn between the spacecraft and Sun, 
limb-darkening, and foreshortening.  

## Velocity Corrections

Dopplergrams are maps of solar surface velocity that must be corrected for the relative
spacecraft velocity, and the velocity due to solar rotation.  

**Spacecraft Velocity Correction:**   

$v_{sc, ij} = - \frac{\delta w_{ij} v_{sc, w_{ij}} + \delta n_{ij} v_{sc, n{ij}} + \delta r_{ij} v_{sc, r_{ij}}}{d_{ij}}$

where $d_{ij} = \sqrt{\delta w_{ij} ^ 2 + \delta n_{ij} ^ 2 + \delta r_{ij} ^ 2 }$.  

```python
def spacecraft_vel(deltaw, deltan, deltar, dij, vmap):
    """
    function to calculate pixel-wise spacecraft velocities for Sunpy map

    Parameters
    ----------
    deltaw: relative westward position of pixel
    deltan: relative northward position of pixel
    deltar: relative radial position of pixel
    dij: distance between pixel ij and spacecraft
    vmap: Sunpy map object (Dopplergram)

    Returns
    -------
    vsc: array of spacecraft velocities

   """

    # velocity of spacecraft relative to sun
    vscw = vmap.meta['obs_vw']
    vscn = vmap.meta['obs_vn']
    vscr = vmap.meta['obs_vr']

    # pixel-wise magnitude of spacecraft velocity
    vsc = - (deltaw * vscw + deltan * vscn + deltar * vscr) / dij

    return vsc
```

**Solar Rotational Velocity:**  

We calculate the solar differential rotation profile based on Snodgrass & Ulrich (1990).  

$\omega(\theta) = \alpha_1 + \alpha_2 sin^2{\theta} + \alpha_3 sin^4{\theta}$   

$\theta$ is latitude and the $\alpha$ values are $14.713, -2.396, -1.787$ respectively.  

```python
def solar_rot_vel(wij, nij, rij, deltaw, deltan, deltar, dij, vmap, a_parameters=[14.713, -2.396, -1.787]):
    """
    function to calculate pixel-wise velocities due to solar rotation

    Parameters
    ----------
    wij: array of westward values for image
    nij: array of northward values for image
    rij: array of radius values for image
    deltaw: relative westward position of pixel
    deltan: relative northward position of pixel
    deltar: relative radial position of pixel
    dij: distance between pixel ij and spacecraft
    vmap: Sunpy map object (Dopplergram)
    a_parameters: array of solar differential rotation parameters from Snodgrass & Ulrich (1990).

    Returns
    -------
    vrot: array of solar rotation velocities

    """

    # apply to cartesian coordinates
    x1 = wij
    y1 = nij * np.cos(np.deg2rad(vmap.meta['crlt_obs'])) + rij * np.sin(np.deg2rad(vmap.meta['crlt_obs']))
    z1 = - nij * np.sin(np.deg2rad(vmap.meta['crlt_obs'])) + rij * np.cos(np.deg2rad(vmap.meta['crlt_obs']))

    hx = x1 * np.cos(np.deg2rad(vmap.meta['crln_obs'])) + z1 * np.sin(np.deg2rad(vmap.meta['crln_obs']))
    hy = y1
    hz = -x1 * np.sin(np.deg2rad(vmap.meta['crln_obs'])) + z1 * np.cos(np.deg2rad(vmap.meta['crln_obs']))

    # apply parameters to determine vrot for given image pixel
    w = (a_parameters[0] + a_parameters[1] * ((np.sin(hy)) ** 2) + a_parameters[2] * (
            (np.sin(hy)) ** 4)) * 1. / 86400. * np.pi / 180.

    # get projection of solar rotation
    vx_rot = w * hz * vmap.meta['rsun_ref']
    vy_rot = 0.
    vz_rot = -w * hx * vmap.meta['rsun_ref']

    v1 = np.cos(np.deg2rad(vmap.meta['crln_obs'])) * vx_rot - np.sin(np.deg2rad(vmap.meta['crln_obs'])) * vz_rot
    v2 = vy_rot
    v3 = np.sin(np.deg2rad(vmap.meta['crln_obs'])) * vx_rot + np.cos(np.deg2rad(vmap.meta['crln_obs'])) * vz_rot

    # project into correct direction
    vrotw = v1
    vrotn = v2 * np.cos(np.deg2rad(vmap.meta['crlt_obs'])) - v3 * np.sin(np.deg2rad(vmap.meta['crlt_obs']))
    vrotr = v2 * np.sin(np.deg2rad(vmap.meta['crlt_obs'])) + v3 * np.cos(np.deg2rad(vmap.meta['crlt_obs']))

    # get full rotational velocity
    vrot = (deltaw * vrotw + deltan * vrotn + deltar * vrotr) / dij

    return vrot
```

## Intensity Corrections

Intensitygrams are continuum filtergrams that need to be corrected for limb-darkening. We 
account for limb-darkening and flatten the continuum using the limb-darkening polynomial
calculated in Allen's Astrophysical Quantities (1993).  

$L_{ij} = 1 - u_2 - v_2 + u_2 cos \theta + v_2 cos^2 \theta$  

The functions for limb-darkening corrections can be found [here](https://github.com/shalverson/NEID_Solar_analysis/blob/master/tamar/tools/lbc_funcs.py).  

## Magnetic Field Corrections

Magnetograms (line-of-sight longitudinal magnetic field maps) are used to determine magnetic 
observables that show strong correlation with RV variation and velocity components. Due to foreshortening, 
the observed magnetic field (magnetogram data) is less than the true radial field by a factor of $\mu = cos(\theta)$. 

We correct for this as seen here:   
$B_{r, ij} = B{obs, ij} / \mu_{ij}$. 

Additionally, we set all magnetic pixels with magnitude lower than the noise level (8G) to 0. 

