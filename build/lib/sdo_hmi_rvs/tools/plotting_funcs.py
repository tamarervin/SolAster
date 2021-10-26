"""
Tamar Ervin
Date: June 23, 2021

functions for plotting images and various
solar observable comparison plots

"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy import stats
from sdo_hmi_rvs.tools.settings_tamar import PlotDir


def scatter_hist(x, y, bins):
    """
    function to create a scatter plot with associated histograms

    Parameters
    ----------
    x: 2D array of x values
    y: 2D array of y values
    bins: number of bins (integer)

    Returns
    -------
    scatter_hist: scatter plot with corresponding histograms

    """

    # make sure bins value is an integer
    if type(bins) != int:
        bins = int(bins)

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.005

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]

    # setup the plot
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_axes(rect_scatter)
    ax_histx = fig.add_axes(rect_histx, sharex=ax)
    ax_histy = fig.add_axes(rect_histy, sharey=ax)

    # no labels
    ax_histx.tick_params(axis="x", labelbottom=False)
    ax_histy.tick_params(axis="y", labelleft=False)

    # scatter plot
    ax.scatter(x, y, edgecolors='k')

    # plot histograms
    ax_histx.hist(x, bins=bins)
    ax_histy.hist(y, bins=bins, orientation='horizontal')


def animate(i, fig, axs, aia):
    """
    function to animate frames to get a movie

    Parameters
    ----------
    i: index
    fig: figure number
    axs: axes Image object
    aia: image to plot

    Returns
    -------
    plot_obj: array of plot objects for FuncAnimation function

    """
    fig.suptitle("AIA %s" % (str(aia[0][i].meta['t_obs'])))
    axs[0, 0].set_title("AIA %s" % (str(aia[0][i].meta['wavelnth'])))
    axs[0, 1].set_title("AIA %s" % (str(aia[1][i].meta['wavelnth'])))
    axs[1, 0].set_title("AIA %s" % (str(aia[2][i].meta['wavelnth'])))
    axs[1, 1].set_title("AIA %s" % (str(aia[3][i].meta['wavelnth'])))
    plot_obj0 = aia[0][i].plot(axes=axs[0, 0])
    # plot_obj0.set_data(aia[0][i].data)
    plot_obj1 = aia[1][i].plot(axes=axs[0, 1])
    # plot_obj1.set_data(aia[1][i].data)
    plot_obj2 = aia[2][i].plot(axes=axs[1, 0])
    # plot_obj2.set_data(aia[2][i].data)
    plot_obj3 = aia[3][i].plot(axes=axs[1, 1])
    # plot_obj3.set_data(aia[3][i].data)
    return plot_obj0, plot_obj1, plot_obj2, plot_obj3


def plot_image(los_image, nfig=None, cmap='gray', title=None):
    """
    function to plot AIA los disk images

    Parameters
    ----------
    los_image: 2D line of sight image to plot
    nfig: figure number
    cmap: colormap to use
    title: figure title

    Returns
    -------

    """

    plot_arr = los_image.data

    plot_arr[plot_arr < .001] = .001

    norm_max = max(1.01, np.nanmax(plot_arr))
    norm = colors.LogNorm(vmin=1.0, vmax=norm_max)

    # plot the initial image
    if nfig is None:
        cur_figs = plt.get_fignums()
        if not nfig:
            nfig = 0
        else:
            nfig = cur_figs.max() + 1

    plt.figure(nfig)

    plt.imshow(plot_arr, extent=[los_image.x.min(), los_image.x.max(), los_image.y.min(), los_image.y.max()],
               origin="lower", cmap=cmap, aspect="equal", norm=norm)
    plt.xlabel("Latitude")
    plt.ylabel("Carrinton Longitude")
    if title is not None:
        plt.title(title)

    return None


def plot_timeseries(x, y, title=None, xlabel=None, ylabel=None, save_fig=None):
    """
    simple mpl function to plot relevant observables

    Parameters
    ----------
    x: x data to plot
    y: y data to plot
    xlabel: label for x axes
    ylabel: label for y axes
    title: title of plot
    save_fig: path to save file, None if not saving

    Returns
    -------

    """

    # plot style
    plot_style = os.path.join(PlotDir.MPL, 'timeseries.mplstyle')
    plt.style.use(plot_style)

    # plot data
    plt.rcParams['figure.figsize'] = [12, 4]
    plt.scatter(x, y, color='cornflowerblue', edgecolors='k', linewidths=1.0)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # save if needed
    if save_fig is not None:
        plt.savefig(save_fig)


def vert_comp_timeseries(x, y_list, title, xlabel=None, ylabel_list=None, save_fig=None):
    """
    function to plot stacked timeseries for comparison

    Parameters
    ----------
    x: x data to plot
    y_list: list of y data to plot
    title: title of plot
    xlabel: label for x axes
    ylabel_list: list of labels for y axes
    save_fig: path to save file, None if not saving

    Returns
    -------

    """

    # plot style
    plot_style = os.path.join(PlotDir.MPL, 'timeseries.mplstyle')
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.cal'] = 'Helvetica'
    plt.style.use(plot_style)

    # set up figure
    fig, axs = plt.subplots(len(y_list), 1, sharex='all', figsize=[12, 3*len(y_list)], gridspec_kw={'hspace': 0.0})
    if title is not None:
        fig.suptitle(title)

    # set up axes labels
    for i in range(0, len(axs)):
        if i == len(axs) - 1:
            axs[i].set(xlabel=xlabel, ylabel=ylabel_list[i])
        else:
            axs[i].set(ylabel=ylabel_list[i])

    # plot data
    color_list = ['lightcoral', 'lightblue', 'plum', 'pink', 'navajowhite']
    for i in range(0, len(axs)):
        axs[i].scatter(x, y_list[i], color=color_list[i], edgecolors='k', linewidths=1.0)
        # axs[i].locator_params(axis="y", nbins=5)
    # axs[1].scatter(x, y_list[1], color='lightblue', edgecolors='k', linewidths=1.0)

    # set axes tick marks
    axs[i].tick_params(labelbottom=True)

    # align y axis labels
    fig.align_ylabels(axs)

    # save if needed
    if save_fig is not None:
        plt.savefig(save_fig)


def vert_plot_timeseries(x, y_list, title, xlabel=None, ylabel_list=None, save_fig=None):
    """
    function to plot stacked timeseries for comparison

    Parameters
    ----------
    x: x data to plot
    y_list: list of y data to plot
    title: title of plot
    xlabel: label for x axes
    ylabel_list: list of labels for y axes
    save_fig: path to save file, None if not saving

    Returns
    -------

    """

    # plot style
    plot_style = os.path.join(PlotDir.MPL, 'timeseries.mplstyle')
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.cal'] = 'Helvetica Neue LT Pro'
    plt.style.use(plot_style)

    # set up figure
    fig, axs = plt.subplots(len(y_list), 1, sharex='all', figsize=[12, 2 * len(y_list)], gridspec_kw={'hspace': 0.0})
    if title is not None:
        fig.suptitle(title)

    # set up axes labels
    for i in range(0, len(axs)):
        if i == len(axs) - 1:
            axs[i].set(xlabel=xlabel, ylabel=ylabel_list[i])
        else:
            axs[i].set(ylabel=ylabel_list[i])

    # plot data
    color_list = ['lightcoral', 'lightblue', 'plum', 'pink']
    for i in range(0, len(axs)):
        axs[i].plot(x, y_list[i], color=color_list[i], linewidth=1.2)
    # axs[1].plot(x, y_list[1], color='steelblue', linewidth=0.8)

    # set axes tick marks
    axs[i].tick_params(labelbottom=True)

    # align y axis labels
    fig.align_ylabels(axs)

    # save if needed
    if save_fig is not None:
        plt.savefig(save_fig)


def correlation_plot(x, y, title, xlabel, ylabel, save_fig=None):
    """
    function to create correlation plot and calculate correlation coefficient

    Parameters
    ----------
    x: x data to plot
    y: y data to plot
    title: title of plot
    xlabel: label for x axes
    ylabel: label for y axes
    save_fig: path to save file, None if not saving

    Returns
    -------

    """

    # plot style
    plot_style = os.path.join(PlotDir.MPL, 'mystyle.mplstyle')
    plt.style.use(plot_style)

    # calculate correlation coefficient
    fig = plt.figure()
    plt.rcParams['figure.figsize'] = [12, 6]
    ax = fig.add_subplot()
    fig.subplots_adjust(top=0.85)

    correlation = stats.spearmanr(x, y)

    # label correlation coefficients
    ax.text(0.1, 0.95, str(np.around(correlation[0], 2)), horizontalalignment='left',
            verticalalignment='top', transform=ax.transAxes, weight='normal')

    # plot data
    plt.scatter(x, y, color='deepskyblue', edgecolors='k', linewidths=1.0)

    # title and axes
    fig.suptitle(title)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # save if needed
    if save_fig is not None:
        plt.savefig(save_fig)


def two_correlation_plots(x_list, y, title, xlabel_list=None, ylabel=None, save_fig=None):
    """
    function for multiple correlation plots

    Parameters
    ----------
    x_list: list of x data to plot
    y: y data to plot
    title: title of plot
    xlabel_list: list of labels for x axes
    ylabel: label for y axes
    save_fig: path to save file, None if not saving

    Returns
    -------

    """

    # plot style
    plot_style = os.path.join(PlotDir.MPL, 'mystyle.mplstyle')
    plt.style.use(plot_style)

    # calculate correlation coefficient
    corr_one = stats.spearmanr(x_list[0], y)
    corr_two = stats.spearmanr(x_list[1], y)

    # setup the plot
    fig, axs = plt.subplots(1, 2, sharey='all', figsize=[12, 6], gridspec_kw={'wspace': 0.1})
    fig.suptitle(title)
    axs[0].set(xlabel=xlabel_list[0], ylabel=ylabel)
    axs[1].set(xlabel=xlabel_list[1])

    # plot data
    axs[0].scatter(x_list[0], y, color='violet', edgecolors='k', linewidths=1.0)
    axs[1].scatter(x_list[1], y, color='violet', edgecolors='k', linewidths=1.0)

    # label correlation coefficients
    axs[0].text(0.1, 0.95, str(np.around(corr_one[0], 2)), horizontalalignment='left',
                verticalalignment='top', transform=axs[0].transAxes, weight='normal')
    axs[1].text(0.1, 0.95, str(np.around(corr_two[0], 2)), horizontalalignment='left',
                verticalalignment='top', transform=axs[1].transAxes, weight='normal')

    # remove tick marks
    for ax in axs:
        ax.tick_params(labelbottom=True)

    # save if needed
    if save_fig is not None:
        plt.savefig(save_fig)


def four_correlation_plots(x_list, y_list, title, xlabel_list=None, ylabel_list=None, save_fig=None):
    """
    function for multiple correlation plots

    Parameters
    ----------
    x_list: list of x data to plot
    y_list: list of y data to plot
    title: title of plot
    xlabel_list: list of labels for x axes
    ylabel_list: list of labels for y axes
    save_fig: path to save file, None if not saving

    Returns
    -------

    """

    # plot style
    plot_style = os.path.join(PlotDir.MPL, 'mystyle.mplstyle')
    plt.style.use(plot_style)

    # calculate correlation coefficient
    corr_one = stats.spearmanr(x_list[0], y_list[0])
    corr_two = stats.spearmanr(x_list[1], y_list[0])
    corr_three = stats.spearmanr(x_list[0], y_list[1])
    corr_four = stats.spearmanr(x_list[1], y_list[1])

    # setup the plot
    fig, axs = plt.subplots(2, 2, sharey='row', sharex='col', figsize=[12, 6], gridspec_kw={'hspace': 0.1, 'wspace': 0.1})
    fig.suptitle(title)
    axs[0, 0].set(ylabel=ylabel_list[0])
    # axs[0, 1].set(xlabel=xlabel_list[1])
    axs[1, 0].set(xlabel=xlabel_list[0], ylabel=ylabel_list[1])
    axs[1, 1].set(xlabel=xlabel_list[1])

    # plot data
    axs[0, 0].scatter(x_list[0], y_list[0], color='violet', edgecolors='k', linewidths=1.0)
    axs[0, 1].scatter(x_list[1], y_list[0], color='violet', edgecolors='k', linewidths=1.0)
    axs[1, 0].scatter(x_list[0], y_list[1], color='violet', edgecolors='k', linewidths=1.0)
    axs[1, 1].scatter(x_list[1], y_list[1], color='violet', edgecolors='k', linewidths=1.0)

    # label correlation coefficients
    axs[0, 0].text(0.1, 0.95, str(np.around(corr_one[0], 2)), horizontalalignment='left',
                verticalalignment='top', transform=axs[0, 0].transAxes, weight='normal')
    axs[0, 1].text(0.1, 0.95, str(np.around(corr_two[0], 2)), horizontalalignment='left',
                verticalalignment='top', transform=axs[0, 1].transAxes, weight='normal')
    axs[1, 0].text(0.1, 0.95, str(np.around(corr_three[0], 2)), horizontalalignment='left',
                verticalalignment='top', transform=axs[1, 0].transAxes, weight='normal')
    axs[1, 1].text(0.1, 0.95, str(np.around(corr_four[0], 2)), horizontalalignment='left',
                verticalalignment='top', transform=axs[1, 1].transAxes, weight='normal')

    # remove tick marks
    # for ax in axs:
    #     ax.tick_params(labelbottom=True)

    # save if needed
    if save_fig is not None:
        plt.savefig(save_fig)


def scatter_overlaid_timeseries(x, y_list, title, xlabel, ylabel_list, save_fig=None):
    """
    function to overlay scatter plots

    Parameters
    ----------
    x: x data to plot
    y_list: list of y data to plot
    title: title of plot
    xlabel: label for x axes
    ylabel_list: list of labels for y axes
    save_fig: path to save file, None if not saving

    Returns
    -------

    """

    # plot style
    plot_style = os.path.join(PlotDir.MPL, 'mystyle.mplstyle')
    plt.style.use(plot_style)

    # plot data
    plt.rcParams['figure.figsize'] = [12, 6]

    fig, ax1 = plt.subplots()

    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel_list[0])
    ax1.set_title(title)
    ax1.scatter(x, y_list[0], color='forestgreen', label=ylabel_list[0])
    ax1.tick_params(axis='y')
    plt.legend()

    ax2 = ax1.twinx()

    ax2.set_ylabel(ylabel_list[1])
    ax2.scatter(x, y_list[1], color='darkviolet', label=ylabel_list[1])
    ax2.tick_params(axis='y')

    plt.title(title)
    plt.legend()

    # save if needed
    if save_fig is not None:
        plt.savefig(save_fig)


