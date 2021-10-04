# Calculation of Solar Magnetic Observables

Using HMI magnetograms, we are able to calculate magnetic observables that
both correlate with each other as well as the velocity features. This shows 
how the magnetic activity is driving the radial velocity variations.  

## Filling Factor

**Disk averaged filling factors of sunspots and plage which gives
the percentage of magnetically active pixels.**  

The filling factor is a calculation of the percentage of active pixels on the solar
surface at any one time. We calculate three different factors: $f_{spot}$, $f_{bright}$,
and $f$. These are the filling factors due to sunspots, faculae/network, and the total
filling factor which is the sum of the two.  

We return the filling factor as a percentage:  
$f = \frac{1}{N_{pix}} \sum_{ij} W_{ij} * 100$

```python
def filling_factor(mu, mmap, active, fac_inds, spot_inds, mu_cutoff=0.3):
    """
    function to calculate filling factor

    Parameters
    ----------
    mu: array of mu (cosine theta) values
    mmap: UNCORRECTED Sunpy map object (Magnetogram)
    active: weights array where active pixels have weight = 1
    fac_inds: array of indices where faculae are detected
    spot_inds: array of indices where sunspots are detected
    mu_cutoff: minimum mu cutoff value

    Returns
    -------
    f_bright: filling factor (%) for bright areas (faculae)
    f_spot: filling factor (%) for dark areas (sunspots)
    f_total: filling factor (%) for timestamp

    """

    # get good mu values
    good_mu = np.where(mu > mu_cutoff)

    # get number of pixels
    npix = len(mmap.data[good_mu])

    # faculae
    faculae = np.zeros(mmap.data.shape)
    faculae[fac_inds] = 1.
    f_bright = np.sum(faculae) / npix * 100

    # sunspots
    spots = np.zeros(mmap.data.shape)
    spots[spot_inds] = 1.
    f_spot = np.sum(spots) / npix * 100

    # get filling factor
    f_total = np.sum(active) / npix * 100

    return f_bright, f_spot, f_total
```

## Unsigned Magnetic Flux

**The disc-averaged, line-of-sight unsigned (unpolarized) magnetic flux of the Sun.**

$|\hat{B_{obs}}| = \frac{\sum_{ij} |B_{obs, ij}| I_{ij}} {\sum_{ij} I_{ij}}$

```python
def unsigned_flux(map_mag_obs, imap):
    """
    calculate unsigned magnetic flux

    Parameters
    ----------
    map_mag_obs: corrected observed magnetic field strength Sunpy map object (Magnetogram)
    imap: UNCORRECTED Sunpy map object (Intensitygram)

    Returns
    -------
    unsign_flux: unsigned magnetic flux

    """

    unsign_flux = np.nansum(np.abs(map_mag_obs.data) * imap.data) / np.nansum(imap.data)

    return np.abs(unsign_flux)
```


