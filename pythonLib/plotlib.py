# Anthony Ho, ahho@stanford.edu, 1/15/2015
# Last update 2/27/2015
"""Python module containing some handy plotting tools"""


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


# Handy function to set the commonly used plot modifiers and apply to plot/figure
def setproperties(fig=None, ax=None, figsize=(10, 10), title=None,
                  legend=None, legendloc=1, legendwidth=2.5,
                  xlabel=None, ylabel=None, xlim=None, ylim=None,
                  scix=False, sciy=False, scilimitsx=(-3, 3), scilimitsy=(-3, 3),
                  logx=False, logy=False, borderwidth=2.5, fontsize=None,
                  legendfontsize=20, tickfontsize=20, labelfontsize=20, titlefontsize=18,
                  xticklabelrot=None, yticklabelrot=None):
    """ Convenient tool to set properties of a plot in a single command"""
    # Get figure and axis handles
    if not fig:
        fig = plt.gcf()
    if not ax:
        ax = plt.gca()

    # Set background color to white
    fig.patch.set_facecolor('w')

    # Define figure size
    plt.rcParams['figure.figsize'] = figsize

    # Set title if provided
    if title:
        ax.set_title(title, y=1.02)
    # Show legend if requested
    if legend:
        legend = plt.legend(loc=legendloc, numpoints=1, fontsize=legendfontsize)
        legend.get_frame().set_linewidth(legendwidth)
    # Set x and y labels if provided
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    # Set x and y limits if provided
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)
    # Apply scientific notation to x and y tick marks if requested
    if scix:
        ax.ticklabel_format(axis='x', style='sci', scilimits=scilimitsx)
    if sciy:
        ax.ticklabel_format(axis='y', style='sci', scilimits=scilimitsy)
    # Change axis to log scale if requested
    if logx:
        ax.set_xscale('log')
    if logy:
        ax.set_yscale('log')

    # Rotate x and y tick labels
    if xticklabelrot:
        xticklabels = ax.get_xticklabels() 
        for ticklabel in xticklabels:
            ticklabel.set_rotation(xticklabelrot)
    if yticklabelrot:
        yticklabels = ax.get_yticklabels() 
        for ticklabel in yticklabels:
            ticklabel.set_rotation(yticklabelrot)

    # Set borderwidth
    plt.rc('axes', linewidth=borderwidth)

    # Set individual fontsizes if fontsize is not specified
    if fontsize is None:
        plt.rc('xtick', labelsize=tickfontsize)
        plt.rc('ytick', labelsize=tickfontsize)
        ax.xaxis.label.set_fontsize(labelfontsize)
        ax.yaxis.label.set_fontsize(labelfontsize)
        ax.title.set_fontsize(titlefontsize)
    # Set all fontsizes to fontsize if fontsize is specified
    else:
        plt.rc('xtick', labelsize=fontsize)
        plt.rc('ytick', labelsize=fontsize)
        ax.xaxis.label.set_fontsize(fontsize)
        ax.yaxis.label.set_fontsize(fontsize)
        ax.title.set_fontsize(fontsize)


# Handy function to make a plot look better
# Deprecated - use setproperties() instead
def makepretty(ax):
    """Make a plot look better"""
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


# Plot a rank-sorted plot from a sequence
def rankSortedPlot(series, ascending=False, name=None,
                   markersize=6, markeredgecolor='b',
                   _useCurrFig=False, **kwargs):
    """Plot rank-sorted plot from a sequence"""
    # Make into a Pandas series if not already one
    if not isinstance(series, pd.core.series.Series):
        series = pd.Series(series)
    # Get name from Pandas series if name is not provided
    if series.name and not name:
        name = series.name

    # Plot the rank-sorted data
    sortedSeries = series.order(ascending=ascending)
    plt.plot(series.index, sortedSeries, marker='.', linestyle='None',
             markersize=markersize, markeredgecolor=markeredgecolor, label=name)
    ax = plt.gca()

    # Set properties
    setproperties(**kwargs)

    if not _useCurrFig:
        plt.show(block=False)

    return ax


# Plot multiple rank-sorted plots in the same plot
# Wrapper of rankSortedPlot
def rankSortedPlotMultiple(listSeries, listName=None, colormap='gist_rainbow',
                           legendloc=1, legendfontsize=20, legendwidth=2.5, **kwargs):
    """Plot rank-sorted plots of multiple sequences"""
    # Initialize
    kwargs['_useCurrFig'] = True
    kwargs.pop('markeredgecolor', None)
    kwargs.pop('name', None)
    nSeries = len(listSeries)
    if not listName:
        listName = [None] * nSeries

    # Define color map
    cm = plt.get_cmap(colormap)

    # Go through the list of sequences
    for i in range(len(listSeries)):
        rankSortedPlot(listSeries[i], name=listName[i],
                       markeredgecolor=cm(1.*i/nSeries), **kwargs)

    # Show legend if requested
    if any(listName):
        legend = plt.legend(loc=legendloc, numpoints=1, fontsize=legendfontsize)
        legend.get_frame().set_linewidth(legendwidth)

    plt.show(block=False)

    return


# Plot box and/or violin plots of variant activities along with their counts
def plotActCount(listVarDF, field, listName=None,
                 transform=None, ref=None, mode='b', bootstrap=2000,
                 logAct=False, logCount=True,
                 actLabel=None, actLim=None,
                 barwidth=0.4, xticklabelrot=None,
                 figsize=None, figunitheight=4, figunitwidth=2, 
                 _xticks=None, **kwargs):
    """Plot box and/or violin plots of variant activities along with their counts"""
    # Define parameters
    numVar = len(listVarDF)
    numSubplots = np.sum(['b' in mode, 'v' in mode]) + 1
    if _xticks is None:
        _xticks = np.arange(1, numVar+1)
    xlim = (0.5, numVar+0.5)

    # Get variants' names from annotation if listName is not provided
    if not listName:
        listName = [varDF['annotation'].iloc[0] for varDF in listVarDF]  ## varDF needs to have its index reset prior
    
    # Transform the variants' activities according to the transform function 
    # if it is provided
    if transform:
        listAct = [transform(varDF[field].values) for varDF in listVarDF]
    else:
        listAct = [varDF[field].values for varDF in listVarDF]

    # Get the variants' count
    listCount = [varDF['count'].iloc[0] for varDF in listVarDF]

    # Make figure with subplots
    if not figsize:
        figwidth = figunitwidth * numVar
        figheight = figunitheight * numSubplots
        figsize = (figwidth, figheight)
    fig, axes = plt.subplots(numSubplots, 1, sharex=True, figsize=figsize)
    fig.patch.set_facecolor('w')
    axCount = axes[-1]
    if 'b' in mode and 'v' in mode:
        axBox = axes[0]
        axViolin = axes[1]
    elif 'b' in mode and not 'v' in mode:
        axBox = axes[0]
    elif 'b' not in mode and 'v' in mode:
        axViolin = axes[0]

    # Make box plot
    if 'b' in mode:
        if ref:
            axBox.plot(range(0,10), [0.032]*10, linestyle='--', color='k')
        sns.boxplot(listAct, positions=_xticks, labels=listName,
                    bootstrap=bootstrap, showmeans=True, widths=barwidth, ax=axBox)
        setproperties(ax=axBox, logy=logAct, ylabel=actLabel, ylim=actLim, **kwargs)

    # Make violin plot
    if 'v' in mode:
        if ref:
            axViolin.plot(range(0,10), [0.032]*10, linestyle='--', color='k')
        sns.violinplot(listAct, positions=_xticks, widths=barwidth, alpha=0.8, ax=axViolin)
        setproperties(ax=axViolin, logy=logAct, ylabel=actLabel, ylim=actLim, **kwargs)

    # Make count bar chart
    axCount.bar(_xticks-barwidth/2, listCount, width=barwidth)
    axCount.set_xticklabels(listName)
    setproperties(ax=axCount, logy=logCount, ylabel='Number of clusters',
                  ylim=(1, None), xlim=xlim, xticklabelrot=xticklabelrot, **kwargs)

    plt.show(block=False)


# Plot the activities and counts of variants given a df of all clusters and list of annotations
# Wrapper of plotActCount
def plotVariants(df, listAnnt, field='params2', transform=None, unit='min', **kwargs):
    """Plot the activities and counts of variants given a df of all clusters and list of annotations
    This function defaults to plotting time/rate constants, but can be used to plot something else too
    """
    # Define unit. Default is min
    if unit in ['second', 's', 'sec']:
        unitTime = 1
    else:
        unitTime = 60

    # Define some commonly used transform functions
    if transform == 'kobs':
        def transformFunc(x): return 1. / x * unitTime
        kwargs['logAct'] = False
        kwargs['actLabel'] = r'$\mathrm{\mathsf{k_{obs}\ (min^{-1})}}$'
    elif transform == 'logkobs':
        def transformFunc(x): return 1. / x * unitTime
        kwargs['logAct'] = True
        kwargs['actLabel'] = r'$\mathrm{\mathsf{k_{obs}\ (min^{-1})}}$'
    elif transform == 'kobsfold':
        def transformFunc(x): return kwargs['ref'] / (1. / x * unitTime)
        kwargs['logAct'] = False
        kwargs['actLabel'] = r'$\mathrm{\mathsf{k_{obs}\ fold change}}$'
    elif transform == 'logkobsfold':
        def transformFunc(x): return kwargs['ref'] / (1. / x * unitTime)
        kwargs['logAct'] = True
        kwargs['actLabel'] = r'$\mathrm{\mathsf{k_{obs}\ fold change}}$'
    else:
        transformFunc = transform
    
    # Get the dataframe for each name in listName
    if df.index.name == 'annotation':
        listVarDf = [df.loc[annt].reset_index() for annt in listAnnt]
    else:
        df2 = df.set_index('annotation')
        listVarDf = [df2.loc[annt].reset_index() for annt in listAnnt]

    # Call plotActCount to plot
    plotActCount(listVarDf, field, transform=transformFunc, **kwargs)


# Plot the activities and counts of single mutants given a df of all clusters and list of positions
# Wrapper of plotVariants, which is a wrapper of pltoSingleMutants
def plotSingleMutants(df, consensus, listPos, field='params2', show='m', **kwargs):
    """Plot the activities and counts of single mutants given a df of all clusters, name of consensus 
    sequence, and list of positions
    This function defaults to plotting time/rate constants, but can be used to plot something else too"""

    
    return 
