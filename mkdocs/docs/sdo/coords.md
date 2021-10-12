# Coordinate Transformations

SDO/HMI Data Products are returned in Helioprojective Cartesian pixel coordinates. In 
order to calculate the spacecraft and solar rotation velocity corrections, we need to
convert this to the Heliographic Carrington Coordinate System (sun centered coordinate system). 
In addition to being sun-centered, Heliographic Carrington coordinates also rotate with
respect to the Carrington rotation period. Transformations are based off Thompson (2006). 

Coordinate functions are [here](https://github.com/tamarervin/sdo_hmi_rvs/blob/main/sdo_hmi_rvs/tools/coord_funcs.py).

## Pixel Array Computation

The first step in the coordinate transformation is computation of the pixel-wise
coordinate arrays with respect to solar center. We do this by taking into account the
reference pixels from the header file. 

We then calculate the array of $\mu$ and radius values for that image, where $\mu = cos(\theta)$. 
First, we determine pixel size in physical coordinates (solar radii) and then calculate
the center-to-limb angles. We get the radius by determining the distance in solar radii 
to each pixel. 

```python
def coordinates(smap):
    """
    calculate array of mu values and cartesian
    coordinates for image

    Parameters
    ----------
    smap: Sunpy map object

    Returns
    -------
    x: array of x coordinates in helioprojective cartesian system
    y: array of y coordinates in helioprojective cartesian system
    pd: array of pixel distances
    r: array of solar radius values
    d: observer distance from sun in solar radii
    mu: array of mu values


    """

    # distance of observer to sun in solar radii
    d = smap.fits_header['DSUN_OBS'] / smap.fits_header['RSUN_REF']

    # focal length in pixels
    f = 180. * 3600. / np.pi / smap.fits_header['CDELT1']

    # get cartesian x, y map scales
    x, y = get_map_scales(smap)

    # distance (in pixels) to pixel
    pd = np.sqrt(x ** 2 + y ** 2)

    # distance (in solar r) to pixel
    pr = f * f * pd * pd + pd ** 4 - d * d * pd ** 4 + 0.J
    r = (d * f * pd - np.sqrt(pr)) / (f * f + pd * pd)

    # separate complex parts
    r = r.real

    # get mu array
    pr = 1 - r ** 2 + 0.J
    cos_alpha = (np.sqrt(pr)).real
    sin_alpha = r.real
    cos_theta = ((d - np.sqrt(pr)) / np.sqrt(r ** 2 + (d - np.sqrt(pr)) ** 2)).real
    sin_theta = (np.sqrt(1 - cos_theta ** 2)).real
    mu = cos_alpha * cos_theta - sin_alpha * sin_theta

    return x, y, pd, r, d, mu
```


## Coordinate Transform

We are able to apply the coordinate transformation based on a calculated rotation 
matrix. We get an array of coordinates in the westward, northward, and radial direction. 

```python
def vel_coords(x, y, pd, r, smap):
    """
    calculate coordinate transformation to heliographic Carrington coordinates

    Parameters
    ----------
    x: array of x coordinates in helioprojective cartesian system
    y: array of y coordinates in helioprojective cartesian system
    pd: array of pixel distances
    r: array of solar radius values
    smap: Sunpy map object

    Returns
    -------
    wij: array of pixel coordinates relative to solar center in westward direction
    nij: array of pixel coordinates relative to solar center in northward direction
    rij: array of pixel coordinates relative to solar center in radial direction

    """

    head = smap.fits_header
    crota2 = head['CROTA2']  # deg

    # transform each pixel to get into Heliographic CR Coordinates
    dw = y * np.sin(np.deg2rad(crota2)) + x * np.cos(np.deg2rad(crota2))
    dn = y * np.cos(np.deg2rad(crota2)) - x * np.sin(np.deg2rad(crota2))

    # get cartesian coordinates for velocity calculations
    pr = 1 - r ** 2 + 0.J
    wij = r * dw / pd
    nij = r * dn / pd
    rij = (np.sqrt(pr)).real

    return wij, nij, rij
```

## Fix Mu Values

Lastly, we remove all pixels with mu values less than 0.1 because our limb-brightening model
does poorly in correcting near the edge of the Solar disk. 


```python
def fix_mu(mu, smaps, mu_cutoff=0.3):
    """
    function to remove pixel values where mu is less than 0.1

    Parameters
    ----------
    mu: cosine of the center to limb angle
    smaps: list of Sunpy map object
    mu_cutoff: minimum mu cutoff value

    Returns
    -------
    mu: corrected cosine of the center to limb angle
    smap: corrected Sunpy map object

    """

    # remove pixel values where mu < mu_cut (0.3)
    bad_inds = np.where(mu <= mu_cutoff)

    # fix arrays
    for smap in smaps:
        smap.data[bad_inds] = 0

    return smaps

```