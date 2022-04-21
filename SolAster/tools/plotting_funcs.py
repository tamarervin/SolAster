"""
Tamar Ervin
Date: June 23, 2021

functions for plotting images and various
solar observable comparison plots

"""

import os
import numpy as np
import matplotlib.pyplot as plt
from SolAster.tools.settings import Inputs


def hmi_plot(int_map, mag_map, vel_map, fac_inds, spot_inds, mu, save_fig=Inputs.save_fig):
    """
    function to plot diagnostic plots showing HMI images and thresholded maps
    Identical to Figure 1 in Ervin et al. (2022) - Accepted.

    Parameters
    ----------
    int_map: Sunpy map
        flattened intensitygram
    mag_map: Sunpy map
        corrected magntogram
    vel_map: Sunpy map
        spacecraft velocity and differential rotation subtracted Dopplergram
    fac_inds: int, array
        array of indices where faculae are detected
    spot_inds: int, array
        array of indices where sunspots are detected
    mu: float, array
        array of mu (cosine theta) values
    save_fig: str
        path to save file, None if not saving

    Returns
    -------

    """

    # make cute plots
    fig, axs = plt.subplots(2, 2, sharey='row', sharex='col', figsize=[12, 12],
                            gridspec_kw={'hspace': 0.1, 'wspace': 0.1})

    # intensity map
    int_data = int_map.data
    int_data = np.where(int_data == 0, np.nan, int_data)
    axs[0, 0].imshow(int_data, cmap=plt.get_cmap('hinodesotintensity'))
    axs[0, 0].set_title("Flattened Continuum Intensity")
    # fig.colorbar(cm.ScalarMappable(cmap=plt.get_cmap('hinodesotintensity')), ax=axs[0, 0])

    # magnetic field map
    mag_data = np.abs(mag_map.data)
    # mag_data = np.where(mag_data == 0, np.nan, mag_data)
    axs[1, 0].imshow(mag_data, cmap=plt.get_cmap('Purples'))
    axs[1, 0].set_title("Unsigned Magnetic Flux Density")
    # fig.colorbar(cm.ScalarMappable(cmap=plt.get_cmap('Purples')), ax=axs[1, 0])

    # Doppler map
    good_mu = np.logical_and(mu > 0.3, mu != np.nan)
    vel_data = np.full(vel_map.data.shape, np.nan)
    vel_data[good_mu] = vel_map.data[good_mu]
    axs[0, 1].imshow(vel_data, cmap=plt.get_cmap('Greys'))
    axs[0, 1].set_title("Line-of-sight Corrected Doppler Velocity")

    # threshold image
    fac = fac_inds.astype(int)
    spot = 1 - spot_inds.astype(int)
    thresh = spot + fac
    axs[1, 1].imshow(thresh, cmap=plt.get_cmap('bwr'))
    axs[1, 1].set_title("Thresholded Map")

    # tick marks
    for i in range(0, 2):
        for j in range(0, 2):
            axs[i, j].set_xticks([])
            axs[i, j].set_yticks([])

    # save if needed
    if save_fig is not None:
        plt.savefig(save_fig)
