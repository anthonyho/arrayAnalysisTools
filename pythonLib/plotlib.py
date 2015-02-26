# Anthony Ho, ahho@stanford.edu, 1/15/2015
# Last update 1/15/2015
"""Python module containing some handy plotting tools"""


import matplotlib.pyplot as plt
import pandas as pd


# Handy function to make a plot look better
def makepretty(ax):

    # Change figure face color to white
    fig = ax.figure
    fig.patch.set_facecolor('white')

    # Eliminate axis ticks on the right and top
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    # Remove the top and right plot frame lines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Set label and title font size
    ax.xaxis.label.set_fontsize(16)
    ax.yaxis.label.set_fontsize(16)
    ax.title.set_fontsize(16)

    return


# Make rank-sorted plot from a sequence
def rankSortedPlot(series, ascending=False,
                   logx=False, logy=False,
                   xlim=None, ylim=None,
                   name=None, legend=True, legendloc=1, legendfontsize=20,
                   xlabel=None, ylabel=None, labelfontsize=20, tickfontsize=20,
                   markersize=6, borderwidth=1.5, _useCurrFig=False):
    """Plot rank-sorted plot from a sequence"""

    #
    if not isinstance(series, pd.core.series.Series):
        series = pd.Series(series)
    if series.name and not name:
        name = series.name

    #
    sortedSeries = series.order(ascending=ascending)
    plt.plot(series.index, sortedSeries, marker='.', linestyle='None',
             markersize=markersize, label=name)
    ax = plt.gca()

    # This should be changed to some standard lib functions
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if legend:
        plt.legend(loc=legendloc, numpoints=1, fontsize=legendfontsize)
    if logy:
        ax.set_yscale('log')
    if logx:
        ax.set_xscale('log')

    #
    plt.xlim(xlim)
    plt.ylim(ylim)

    #
    plt.rc('axes', linewidth=borderwidth)
    plt.rc('xtick', labelsize=tickfontsize)
    plt.rc('ytick', labelsize=tickfontsize)
    ax.xaxis.label.set_fontsize(labelfontsize)
    ax.yaxis.label.set_fontsize(labelfontsize)

    if not _useCurrFig:
        fig = plt.gcf()
        fig.patch.set_facecolor('w')
        plt.show(block=False)

    return ax


# 
def rankSortedPlotMultiple(listSeries, listnames=None, figsize=(7.5, 7.5), **kwargs):
    """Plot rank-sorted plots of multiple sequences"""
    fig = plt.figure(figsize=figsize)

    #
    for i in range(len(listSeries)):
        if listnames:
            rankSortedPlot(listSeries[i], name=listnames[i], _useCurrFig=True, **kwargs)
        else:
            rankSortedPlot(listSeries[i], _useCurrFig=True, **kwargs)

    #
    fig.patch.set_facecolor('w')
    plt.show(block=False)

    return
