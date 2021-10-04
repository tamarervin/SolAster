# Identification of Solar Regions 

In order to calculate the various pertinent velocity components we need
to differentiate between different types of solar regions. We do two identification
cuts both based on thresholding values from Yeo et al. 2013.  

## Magnetic Thresholding

The first threshold is to identify active regions. We do this with a simple 
magnetic threshold value of 24G. The threshold is applied to the unsigned magnetic
field strength (magnetic field corrected for foreshortening). 

Pixels with magnetic field magnitude above the threshold are marked as active and those below the threshold 
are quiet Sun. Isolated active pixels are re-identified as quiet Sun based on area thresholding. 

Isolated pixels are denoted as any pixel without four connected neighbors. The connection
matrix is built using the <code>skimage.regionsprops</code> function.

```python
def mag_thresh(mu, mmap, Br_cutoff=24, mu_cutoff=0.3):
    """
    function to calculate magnetic threshold and differentiate between magnetically active regions and quiet Sun

    Parameters
    ----------
    mu: array of mu (cosine theta) values
    mmap: corrected (unsigned magnetic field) Sunpy map object (Magnetogram)
    Br_cutoff: minimum cutoff value (in Gauss) for thresholding active regions
    mu_cutoff: minimum mu cutoff value for data to ignore

    Returns
    -------
    active: weights array where active pixels are 1
    quiet: weights array where active pixels are 0

    """

    # get active region indices
    active_inds = np.where(np.abs(mmap.data) * mu > Br_cutoff)
    bad_mu = np.where(mu <= mu_cutoff)

    # make active region array
    active = np.zeros(mu.shape)
    active[active_inds] = 1.
    active[bad_mu] = 0.

    # find isolated pixels
    # get area
    y_labeled = label(active, connectivity=2, background=0)
    y_area = [props.area for props in regionprops(y_labeled)]

    # area constraint
    good_area = np.where(np.array(y_area) > 5)
    good_area = good_area[0] + 1
    active_indices = np.isin(y_labeled, good_area)

    # create weights array
    active[~active_indices] = 0

    # get quiet indices
    quiet = 1 - active

    return active, quiet
```

## Intensity Thresholding

Once we detect active regions, we want to identify sunspots (dark) and faculae/network regions (bright).
These regions can be categorized by intensity and therefore an intensity cutoff value
is calculated relative to the overall mean flattened intensity of the quiet Sun:   

$I_{quiet} = \frac{\sum_{ij} I_{flat, ij} W_{ij}}{\sum_{ij} W_{ij}}$ where $W_{ij}$ 
is the quiet Sun weighting array returned from magnetic thresholding.  

The intensity threshold is then $I_{thresh} = 0.89I_{quiet}$.  

```python
def int_thresh(map_int_cor, active, quiet):
    """
    function to do intensity thresholding and differentiate between faculae (bright) and sunspots (dark)

    Parameters
    ----------
    map_int_cor: corrected (limb-darkening) Sunpy map object (Intensitygram)
    active: weights array where active pixels are 1
    quiet: weights array where active pixels are 0

    Returns
    -------
    fac_inds: array of indices where faculae are detected
    spot_inds: array of indices where sunspots are detected

    """
    # flattened intensity data
    Iflat = map_int_cor.data

    # calculate quiet sun intensity
    int_quiet = np.nansum(Iflat * quiet) / np.nansum(quiet)

    # intensity threshold
    int_cutoff = 0.89 * int_quiet

    # get faculae
    fac_inds = np.logical_and((Iflat > int_cutoff), (active > 0.5))

    # get sunspots
    spot_inds = np.logical_and((Iflat <= int_cutoff), (active > 0.5))

    return fac_inds, spot_inds
```


