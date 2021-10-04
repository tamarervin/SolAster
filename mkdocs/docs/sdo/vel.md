# Calculation of RV Components

After the data corrections and identification of solar regions, we are able to calculate
the individual feature components that make up the SDO derived RV variation. 

## Photometric Velocity 
**Estimation of rotational Doppler imbalance due to active (faculae/darkspot) regions**  

Calculation of the photometric velocity due the presence of dark spots/bright faculae causing a Doppler imbalance. Requires
an active region based magnetic weighting array. A weight of 0 denotes quiet-Sun and
1 denotes active pixel.

The equation for this velocity calculation is:
$v_{phot} = \frac{\sum_{ij} (v_{ij} - \delta v_{sc, ij})(I_{ij} - K) W_{ij}}
 {\sum_{ij} I_{ij}}$.

K is a scaling factor based off the limb-brightening correction.
$K = \frac{\sum_{ij} I_{ij} L_{ij} W_{ij}}{\sum_{ij} L_{ij}^2 W_{ij}}$

In this case $I_{ij}$ is the uncorrected (original) intensity data.  

We are also able to use out intensity weighting to calculate photometric velocities due to 
just faculae/network or sunspot regions. 

```python
def v_phot(quiet, active, Lij, vrot, imap, mu, fac_inds, spot_inds, mu_cutoff=0.3):
    """
    function to calculate photometic velocity due to rotational Doppler variation

    Parameters
    ----------
    quiet: weights array where active pixels have weight = 0
    active: weights array where active pixels have weight = 1
    Lij: limb-darkening polynomial function
    vrot: solar rotational velocity
    imap: UNCORRECTED Sunpy map object (Intensitygram)
    mu: array of mu values
    fac_inds: array of indices where faculae are detected
    spot_inds: array of indices where sunspots are detected
    mu_cutoff: minimum mu cutoff value

    Returns
    -------
    v_phot: photospheric velocity perturbation (float)

    """

    # get good mu values
    good_mu = np.where(mu > mu_cutoff)

    # calculate K scaling factor
    K = np.nansum(imap.data * Lij * quiet) / np.sum((Lij[good_mu] ** 2) * quiet[good_mu])

    # calculate photospheric velocity
    v_phot = np.nansum(np.real(vrot) * (imap.data - K * Lij) * active) / np.nansum(imap.data)

    # faculae driven photospheric velocity
    vphot_bright = np.nansum(np.real(vrot) * (imap.data - K * Lij) * fac_inds) / np.nansum(imap.data)

    # sunspots driven photospheric velocity
    vphot_spot = np.nansum(np.real(vrot) * (imap.data - K * Lij) * spot_inds) / np.nansum(imap.data)

    return v_phot, vphot_bright, vphot_spot
```

## Convective Velocity 
**Estimation due to the suppression of the convective blueshift by magnetically active regions**

Calculation of suppression of convective blueshift due to active regions (mainly faculae).

The basic premise for calculation is to calculate the disc-averaged velocity of the Sun
and subtract from that the quiet-Sun velocity.

We then calculate the convective velocity as what is left:
$v_{conv} = v - v_{quiet}$ 

### Quiet Sun Velocity 
**Estimation of average RV of quiet-Sun due to convective motion**  

Calculate the velocity due to the convective motion of the quiet-Sun. Requires
the magnetic weighting array we built in step six, part one. A weight of 0 denotes an
active pixel, and 1 is quiet-Sun.

The equation for calculating this quiet-Sun velocity is as such:
$v_{quiet} = \frac{\sum_{ij} (v_{ij} - \delta v_{sc, ij} - \delta v_{rot, ij}) I_{ij} W_{ij}}
 {\sum_{ij} I_{ij} W_{ij}}$.

In this case $I_{ij}$ is the uncorrected (original) intensity data.  

```python
def v_quiet(map_vel_cor, imap, quiet):
    """
    function to calculate velocity due to convective motion of quiet-Sun

    Parameters
    ----------
    map_vel_cor: corrected (velocities) Sunpy map object (Dopplergram)
    imap: UNCORRECTED Sunpy map object (Intensitygram)
    quiet: weights array where active pixels have weight = 0

    Returns
    -------
    v_quiet: quiet-Sun velocity (float)

    """

    v_quiet = np.nansum(map_vel_cor.data * imap.data * quiet) / np.nansum(
        imap.data * quiet)

    return v_quiet
```

### Full Disc Velocity 
**Total corrected (spacecraft motion, solar rotation) disk-averaged velocity of the Sun**  

The equation for the disc-averaged velocity is:
$v = \frac{\sum_{ij} (v_{ij} - \delta v_{sc, ij} - \delta v_{rot, ij}) I_{ij}}{\sum_{ij} I_{ij}}$.  

```python
def v_disc(map_vel_cor, imap):
    """
    function to calculate disc-averaged velocity of Sun

    Parameters
    ----------
    map_vel_cor: corrected (velocities) Sunpy map object (Dopplergram)
    imap: UNCORRECTED Sunpy map object (Intensitygram)

    Returns
    -------
    v_disc: disc averaged velocity of Sun (float)

    """

    v_disc = np.nansum(map_vel_cor.data * imap.data) / np.nansum(imap.data)

    return v_disc
```
