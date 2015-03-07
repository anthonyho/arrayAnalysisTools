# Anthony Ho, ahho@stanford.edu, 1/15/2015
# Last update 2/27/2015
"""Python module containing some handy plotting tools"""


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import seqlib


# Handy function to set the commonly used plot modifiers and apply to plot/figure
def setproperties(fig=None, ax=None, figsize=(10, 10),
                  suptitle=None, title=None,
                  legend=None, legendloc=1, legendwidth=2.5,
                  xlabel=None, ylabel=None, xlim=None, ylim=None,
                  scix=False, sciy=False, scilimitsx=(-3, 3), scilimitsy=(-3, 3),
                  logx=False, logy=False, majorgrid=None, minorgrid=None,
                  borderwidth=2.5, pad=1.6,
                  fontsize=None, legendfontsize=20, tickfontsize=20, 
                  labelfontsize=20, titlefontsize=18,
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

    # Set titles if provided
    if suptitle:
        if fontsize is None:
            fig.suptitle(suptitle, fontsize=titlefontsize, y=0.99)
        else:
            fig.suptitle(suptitle, fontsize=fontsize, y=0.99)
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
    # Set major and minor grid in plot
    if majorgrid is not None:
        ax.grid(b=majorgrid, which='major')
    if minorgrid is not None:
        ax.grid(b=minorgrid, which='minor')
    
    # Rotate x and y tick labels if requested
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

    # Set tight figure and padding
    fig.tight_layout(pad=pad)


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
                   _useCurrFig=False, show=True, **kwargs):
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
    fig = plt.gcf()

    # Set properties
    setproperties(**kwargs)

    if not _useCurrFig and show:
        plt.show(block=False)

    return fig, fig.axes


# Plot multiple rank-sorted plots in the same plot
# Wrapper of rankSortedPlot
def rankSortedPlotMultiple(listSeries, listName=None, colormap='gist_rainbow',
                           legendloc=1, legendfontsize=20, legendwidth=2.5,
                           show=True, **kwargs):
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
        fig, ax = rankSortedPlot(listSeries[i], name=listName[i],
                                 markeredgecolor=cm(1.*i/nSeries), **kwargs)

    # Show legend if requested
    if any(listName):
        legend = plt.legend(loc=legendloc, numpoints=1, fontsize=legendfontsize)
        legend.get_frame().set_linewidth(legendwidth)

    if show:
        plt.show(block=False)

    return fig, ax


# Plot box and/or violin plots of variant activities along with their counts
# This is the most low level of the plotActCount-related script that is not 
# intended to be used on its own; consider using plotVariants or 
# plotSingleMutants instead
def plotActCount(listVarDF, field, listName=None, color=None,
                 transform=None, ref=None, plotmode='b', bootstrap=2000, inner=None,
                 logAct=False, logCount=True,
                 actLabel=None, actLim=None,
                 barwidth=0.4, xticklabelrot=None,
                 figsize=None, figunitheight=4, figunitwidth=2, maxfigwidth=32,
                 _xticks=None, show=True, **kwargs):
    """Plot box and/or violin plots of variant activities along with their counts

    Required arguments:
      listVarDF - a list of dataframe, each of them belongs to a single variant
      field - the name of the column of the dataframes to be plotted
    """
    # Define parameters
    numVar = len(listVarDF)
    numSubplots = np.sum(['b' in plotmode, 'v' in plotmode]) + 1
    if _xticks is None:
        _xticks = np.arange(1, numVar+1)
        xlim = (0.5, numVar+0.5)
    else:
        xlim = (0.5, max(_xticks)+0.5)
    if not inner:
        maxNumClusters = max([len(varDF) for varDF in listVarDF])
        if maxNumClusters > 500:
            inner = 'box'
        else:
            inner = 'stick'

    # Get variants' names from annotation if listName is not provided
    if not listName:
        listName = [varDF.index[0] if varDF.index.name == 'annotation'
                    else varDF['annotation'].iloc[0]
                    for varDF in listVarDF]

    # Transform the variants' activities according to the transform function 
    # if it is provided
    if transform:
        listAct = [transform(varDF[field].values) for varDF in listVarDF]
    else:
        listAct = [varDF[field].values for varDF in listVarDF]

    # Get the variants' count
    listRawCount = np.array([varDF['count'].iloc[0] for varDF in listVarDF])
    listFilteredCount = np.array([len(varDF) for varDF in listVarDF])

    # Make figure with subplots
    if not figsize:
        figwidth = min(figunitwidth * numVar, maxfigwidth)
        figheight = figunitheight * numSubplots
        figsize = (figwidth, figheight)
    fig, axes = plt.subplots(numSubplots, 1, sharex=True, figsize=figsize)
    fig.patch.set_facecolor('w')
    axCount = axes[-1]
    if 'b' in plotmode and 'v' in plotmode:
        axBox = axes[0]
        axViolin = axes[1]
    elif 'b' in plotmode and not 'v' in plotmode:
        axBox = axes[0]
    elif 'b' not in plotmode and 'v' in plotmode:
        axViolin = axes[0]

    # Make box plot
    if 'b' in plotmode:
        if ref:
            axBox.plot(xlim, [ref]*2, linestyle='--', color='k')
        sns.boxplot(listAct, positions=_xticks, color=color,
                    showmeans=True, bootstrap=bootstrap,
                    widths=barwidth, ax=axBox)
        setproperties(ax=axBox, logy=logAct, ylabel=actLabel,
                      ylim=actLim, majorgrid=True, **kwargs)

    # Make violin plot
    if 'v' in plotmode:
        if ref:
            axViolin.plot(xlim, [ref]*2, linestyle='--', color='k')
        sns.violinplot(listAct, positions=_xticks, color=color,
                       widths=barwidth, alpha=0.8, inner=inner, ax=axViolin)
        setproperties(ax=axViolin, logy=logAct, ylabel=actLabel,
                      ylim=actLim, majorgrid=True, **kwargs)

    # Make count bar chart
    axCount.bar(_xticks-barwidth/2, listFilteredCount, width=barwidth)
    axCount.bar(_xticks-barwidth/2, listRawCount-listFilteredCount, width=barwidth,
                bottom=listFilteredCount, color=sns.xkcd_rgb["tangerine"])
    axCount.set_xticklabels(listName)
    setproperties(ax=axCount, logy=logCount, ylabel='Number of clusters',
                  ylim=(1, None), xlim=xlim, xticklabelrot=xticklabelrot, **kwargs)

    if show:
        plt.show(block=False)
    
    return fig, fig.axes


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

    # Make the df indexed by annotation if not already
    if df.index.name == 'annotation':
        df2 = df.copy()
    else:
        df2 = df.set_index('annotation')

    # See if user has provided an reference annotation; if so, compute
    # reference activity from reference annotation
    refAnnt = kwargs.get('ref', None)
    computeRefFromAnnt = isinstance(refAnnt, basestring)

    # Define some commonly used transform functions
    if transform == 'kobs':
        def transformFunc(x): return 1. / x * unitTime
        kwargs['logAct'] = False
        kwargs['actLabel'] = r'$\mathrm{\mathsf{k_{obs}\ (min^{-1})}}$'
        if computeRefFromAnnt:
            kwargs['ref'] = transformFunc(df2.loc[[refAnnt]][field]).median()
    elif transform == 'logkobs':
        def transformFunc(x): return 1. / x * unitTime
        kwargs['logAct'] = True
        kwargs['actLabel'] = r'$\mathrm{\mathsf{k_{obs}\ (min^{-1})}}$'
        if computeRefFromAnnt:
            kwargs['ref'] = transformFunc(df2.loc[[refAnnt]][field]).median()
    elif transform == 'kobsfold':
        if computeRefFromAnnt:
            ref = (1. / df2.loc[[refAnnt]][field] * unitTime).median()
        def transformFunc(x): return ref / (1. / x * unitTime)
        kwargs['logAct'] = False
        kwargs['actLabel'] = r'$\mathrm{\mathsf{k_{obs}\ fold\ change}}$'
        kwargs.pop('ref', None)
    elif transform == 'logkobsfold':
        if computeRefFromAnnt:
            ref = (1. / df2.loc[[refAnnt]][field] * unitTime).median()
        def transformFunc(x): return ref / (1. / x * unitTime)
        kwargs['logAct'] = True
        kwargs['actLabel'] = r'$\mathrm{\mathsf{k_{obs}\ fold\ change}}$'
        kwargs['ref'] = 1
    else:
        transformFunc = transform
        if computeRefFromAnnt:
            kwargs['ref'] = transformFunc(df2.loc[[refAnnt]][field]).median()

    # Get the dataframe for each name in listName
    listVarDf = [df2.loc[[annt]] for annt in listAnnt]  # note the double square brackets to make sure output is always a DF to prevent n=1 issue

    # Call plotActCount to plot
    return plotActCount(listVarDf, field, transform=transformFunc, **kwargs)


# Plot the activities and counts of single mutants given a df of all clusters and list of positions
# Wrapper of plotVariants, which is a wrapper of pltoSingleMutants
def plotSingleMutants(df, consensus, listPos, muttype='m',
                      colorbymut=True, fullname=False, collapse=False, **kwargs):
    """Plot the activities and counts of single mutants given a df of all clusters, name of consensus 
    sequence, and list of positions
    This function defaults to plotting time/rate constants, but can be used to plot something else too"""
    # Define constants
    allBases = ['A', 'T', 'C', 'G']
    allColors = {'A': sns.xkcd_rgb["windows blue"],
                 'T': sns.xkcd_rgb["amber"],
                 'C': sns.xkcd_rgb["faded green"],
                 'G': sns.xkcd_rgb["reddish orange"],
                 'x': sns.xkcd_rgb["dusty purple"]}
    defaultColor = sns.xkcd_rgb["greyish"]

    # Compute how many bars to be plotted at each position and xtickpos
    numBarDict = {'m': 3, 'i': 3, 'd': 1}
    numBar = np.sum([numBarDict.get(i, 0) for i in muttype])
    gap = int(numBar) / 3
    xtickspos = [range((numBar+gap)*pos+1, (numBar+gap)*pos+numBar+1) for pos in range(0, len(listPos))]
    xtickspos = np.array([item for sublist in xtickspos for item in sublist])

    # Make the list of annotations to plot
    listAnnt = []
    listName = []
    listColor = []
    for pos in listPos:

        # Add mismatches
        if 'm' in muttype:
            # Get base at the current position
            currBase = pos[0]
            # Get the 3 other bases at the current position
            allOtherBases = seqlib.allOtherBases(currBase)
            # Make the list of annotations of the 3 possible mismatches at the current position
            mAnnts = [consensus+':1:0:0:'+pos+base+':::' for base in allOtherBases]
            mNames = [pos+base for base in allOtherBases]
            mColors = [allColors.get(base, defaultColor) for base in allOtherBases]
            # Add to listAnnt and listName
            listAnnt.extend(mAnnts)
            listName.extend(mNames)
            listColor.extend(mColors)
        
        # Add insertions
        if 'i' in muttype:
            # Make the list of annotations of the 4 possible insertions at the current position
            iAnnts = [consensus+':0:1:0::+'+pos[1:]+base+'::' for base in allBases]
            iName = ['+'+pos[1:]+base for base in allBases]
            iColors = [allColors.get(base, defaultColor) for base in allBases]
            # Add to listAnnt and listName
            listAnnt.extend(iAnnts)
            listName.extend(iName)
            listColor.extend(iColors)

        # Add deletions
        if 'd' in muttype:
            # Make annotation for the only possible deletion at the current position
            dAnnt = consensus+':0:0:1:::'+pos+'x:'
            dName = pos+'x'
            dColor = allColors.get('x', defaultColor)
            # Add to listAnnt and listName
            listAnnt.append(dAnnt)
            listName.append(dName)
            listColor.append(dColor)

    # Check if annotations exist in df, if not remove from listAnnt, listName, listColor, and xtickspos
    if df.index.name == 'annotation':
        listValidAnnt = [i for i, annt in enumerate(listAnnt) if annt in df.index]
    else:
        df2 = df.set_index('annotation')
        listValidAnnt = [i for i, annt in enumerate(listAnnt) if annt in df2.index]

    listAnnt = [item for i, item in enumerate(listAnnt) if i in listValidAnnt]
    listName = [item for i, item in enumerate(listName) if i in listValidAnnt]
    listColor = [item for i, item in enumerate(listColor) if i in listValidAnnt]
    xtickspos = np.array([item for i, item in enumerate(xtickspos) if i in listValidAnnt])

    # Call plotVariants to plot
    if not colorbymut:
        listColor = None
    if fullname:
        listName = None
    if collapse:
        xtickspos = None

    return plotVariants(df, listAnnt, listName=listName, color=listColor, _xticks=xtickspos,
                        xticklabelrot=90, **kwargs)
