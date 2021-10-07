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

## Area Based Magnetic Observables

In addition to the unsigned flux and filling factor calculations described above, we calculate
these quantities with an eye for what types of regions are causing variations in magnetic
activity. The calculations for area thresholded regions is outlined below.

### Area Thresholded Filling Factor

```python
def area_filling_factor(active, area, mu, mmap, fac_inds, athresh=20, mu_cutoff=0.3):
    """
    calculate filling factor for regions thresholded by area
    - differentiate between large and small regions
    - differentiate between plage (large) and network (small) bright regions

    Parameters
    ----------
    active: weights array where active pixels have weight = 1
    area: area of each active region weighted by its intensity
    mu: array of mu (cosine theta) values
    mmap: UNCORRECTED Sunpy map object (Magnetogram)
    fac_inds: array of indices where faculae are detected
    athresh: area threshold value between large and small regions (in uHem)
    mu_cutoff: minimum mu cutoff value for usable data

    Returns
    -------
    f_small: filling factor (%) for small magnetically active regions
    f_large: filling factor (%) for large magnetically active regions
    f_network: filling factor (%) for network (small, bright magnetically active) regions
    f_plage: filling factor (%) for plage (large, bright magnetically active) regions
    f_nonconv: filling factor (%) for regions that do not suppress convective blueshift

    """

    # get good mu values
    good_mu = np.where(mu > mu_cutoff)

    # get number of pixels
    npix = len(mmap.data[good_mu])

    # get quiet pixels
    quiet = 1 - active

    # get filling factor for 'small' magnetic features
    small = np.zeros(mmap.data.shape)
    small_inds = np.logical_and(active > 0.5, area < athresh)
    small[small_inds] = 1.
    f_small = np.nansum(small) / npix * 100

    # get filling factor for 'large' magnetic features
    large = np.zeros(mmap.data.shape)
    large_inds = np.logical_and(active > 0.5, area > athresh)
    large[large_inds] = 1.
    f_large = np.nansum(large) / npix * 100

    # get filling factor for network (small, faculae regions)
    network = np.zeros(mmap.data.shape)
    network_inds = np.logical_and(small > 0.5, fac_inds > 0.5)
    network[network_inds] = 1.
    f_network = np.nansum(network) / npix * 100

    # get filling factor for plage (large, faculae regions)
    plage = np.zeros(mmap.data.shape)
    plage_inds = np.logical_and(large > 0.5, fac_inds > 0.5)
    plage[plage_inds] = 1.
    f_plage = np.nansum(plage) / npix * 100

    # get filling factor for small, non-convective regions
    nonconv = np.zeros(mmap.data.shape)
    nonconv_inds = np.logical_and(quiet > 0.5, small > 0.5)
    nonconv[nonconv_inds] = 1.
    f_nonconv = np.nansum(nonconv) / npix * 100

    return f_small, f_large, f_network, f_plage, f_nonconv
```

### Area Thresholded Unsigned Flux 
```python
def area_unsigned_flux(map_mag_obs, imap, area, active, athresh=20):
    """
    calculate the magnetic flux for different regions based on area cut
    and magnetic activitiy

    Parameters
    ----------
    map_mag_obs: corrected observed magnetic field strength Sunpy map object (Magnetogram)
    imap: UNCORRECTED Sunpy map object (Intensitygram)
    area: area of each active region weighted by its intensity
    active: weights array where active pixels have weight = 1

    Returns
    -------
    quiet_flux: magnetic flux of quiet-Sun regions
    ar_flux: magnetic flux of active regions
    conv_flux: magnetic flux of regions that suppress convective blueshift
    pol_flux: magnetic flux of polarized regions
    pol_conv_flux: magnetic flux of polarized regions that suppress the convective blueshift

    """

    # get data arrays
    i_data = imap.data
    m_data = map_mag_obs.data
    mabs_data = np.abs(m_data)
    quiet = 1 - active

    # get large regions array
    large = np.zeros(m_data.shape)
    large_inds = np.logical_and(active > 0.5, area > athresh)
    large[large_inds] = 1.

    # calculate relevant fluxes
    quiet_flux = np.nansum(mabs_data * i_data * quiet) / np.nansum(i_data * quiet)
    ar_flux = np.nansum(mabs_data * i_data * active) / np.nansum(i_data * active)
    conv_flux = np.nansum(mabs_data * i_data * large) / np.nansum(i_data * large)
    pol_flux = np.nansum(m_data * i_data) / np.nansum(i_data)
    pol_conv_flux = np.nansum(m_data * i_data * large) / np.nansum(i_data * large)

    return quiet_flux, ar_flux, conv_flux, pol_flux, pol_conv_flux
```

