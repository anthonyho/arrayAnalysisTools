# Anthony Ho, ahho@stanford.edu, 1/15/2015
# Last update 3/10/2015
"""Python module containing some handy plotting tools"""


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatter
import seaborn as sns
import pandas as pd
import seqlib


# Handy function to set the commonly used plot modifiers and apply to plot/figure
def setproperties(fig=None, ax=None, figsize=None,
                  suptitle=None, title=None,
                  legend=None, legendloc=1, legendwidth=2.5,
                  xlabel=None, ylabel=None, xlim=None, ylim=None,
                  scix=False, sciy=False, scilimitsx=(-3, 3), scilimitsy=(-3, 3),
                  logx=False, logy=False, majorgrid=None, minorgrid=None,
                  borderwidth=2.5, tight=True, pad=1.6,
                  fontsize=None, legendfontsize=20, tickfontsize=20,
                  labelfontsize=20, titlefontsize=18, suptitlefontsize=20,
                  xticklabelrot=None, yticklabelrot=None):
    """ Convenient tool to set properties of a plot in a single command"""
    # Get figure and axis handles
    if not fig:
        fig = plt.gcf()
    if not ax:
        ax = plt.gca()

    # Set background color to white
    fig.patch.set_facecolor('w')

    # Define figure size if provided
    if figsize:
        fig.set_size_inches(figsize, forward=True)

    # Set titles if provided
    if suptitle is not None:
        if fontsize is None:
            fig.suptitle(suptitle, fontsize=suptitlefontsize, y=0.99)
        else:
            fig.suptitle(suptitle, fontsize=fontsize, y=0.99)
    if title is not None:
        ax.set_title(title, y=1.02)
    # Show legend if requested
    if legend:
        legend = plt.legend(loc=legendloc, numpoints=1, fontsize=legendfontsize)
        legend.get_frame().set_linewidth(legendwidth)
    # Set x and y labels if provided
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    # Set x and y limits if provided
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
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
    if xticklabelrot is not None:
        xticklabels = ax.get_xticklabels()
        for ticklabel in xticklabels:
            ticklabel.set_rotation(xticklabelrot)
    if yticklabelrot is not None:
        yticklabels = ax.get_yticklabels()
        for ticklabel in yticklabels:
            ticklabel.set_rotation(yticklabelrot)

    # Set borderwidth (not visible if using seaborn default theme)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(borderwidth)

    # Set individual fontsizes if fontsize is not specified
    if fontsize is None:
        plt.setp(ax.get_xticklabels(), fontsize=tickfontsize)
        plt.setp(ax.get_yticklabels(), fontsize=tickfontsize)
        ax.xaxis.label.set_fontsize(labelfontsize)
        ax.yaxis.label.set_fontsize(labelfontsize)
        ax.title.set_fontsize(titlefontsize)
    # Set all fontsizes to fontsize if fontsize is specified
    else:
        plt.setp(ax.get_xticklabels(), fontsize=fontsize)
        plt.setp(ax.get_yticklabels(), fontsize=fontsize)
        ax.xaxis.label.set_fontsize(fontsize)
        ax.yaxis.label.set_fontsize(fontsize)
        ax.title.set_fontsize(fontsize)

    # Set tight figure and padding
    if tight:
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
    minorgridLogAct = logAct

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
    elif 'b' in plotmode and 'v' not in plotmode:
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
        setproperties(ax=axBox, logy=logAct, ylabel=actLabel, ylim=actLim,
                      majorgrid=True, minorgrid=minorgridLogAct, **kwargs)

    # Make violin plot
    if 'v' in plotmode:
        axViolin2 = axViolin.twinx()  # Make an invisible second y-axis to deal with issues with plotting violin plots in log space
        if ref:
            axViolin.plot(xlim, [ref]*2, linestyle='--', color='k')
        if logAct:
            listLogAct = [np.log10(act) for act in listAct]
            sns.violinplot(listLogAct, positions=_xticks, color=color,
                           widths=barwidth, alpha=0.8, inner=inner, ax=axViolin2)
        else:
            sns.violinplot(listAct, positions=_xticks, color=color,
                           widths=barwidth, alpha=0.8, inner=inner, ax=axViolin)
        # Make the left y-axis to plot in log scale
        setproperties(ax=axViolin, logy=logAct, ylabel=actLabel, ylim=actLim,
                      majorgrid=True, minorgrid=minorgridLogAct, **kwargs)
        # Make the right y-axis to plot the already logged data in linear scale
        setproperties(ax=axViolin2, ylim=np.log10(actLim))
        axViolin2.get_yaxis().set_visible(False)

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
    elif transform is None:
        def transformFunc(x): return x
        if computeRefFromAnnt:
            kwargs['ref'] = transformFunc(df2.loc[[refAnnt]][field]).median()
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
def plotSingleMutants(df, consensus, listSeqPos, muttype='m',
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
    xtickspos = [range((numBar+gap)*pos+1, (numBar+gap)*pos+numBar+1) for pos in range(0, len(listSeqPos))]
    xtickspos = np.array([item for sublist in xtickspos for item in sublist])

    # Make the list of annotations to plot
    listAnnt = []
    listName = []
    listColor = []
    for pos in listSeqPos:

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


# Draw a box around a specific cell in a heatmap
def _drawBox(ax, xpos, ypos, color='g', linewidth=4):
    ax.plot([xpos, xpos], [ypos, ypos+1], color=color, linewidth=linewidth)
    ax.plot([xpos+1, xpos+1], [ypos, ypos+1], color=color, linewidth=linewidth)
    ax.plot([xpos, xpos+1], [ypos, ypos], color=color, linewidth=linewidth)
    ax.plot([xpos, xpos+1], [ypos+1, ypos+1], color=color, linewidth=linewidth)
    return


# High level function to plot base triple heat maps
def plotBaseTriples(dfUnqClusters, dfMut, field='params2.median', transform=None, ref=None,
                    orderBaseTriple=None, figsize=(12, 12), suptitle=None, actLabel=None,
                    vmin=None, vmax=None, cmap=None, c_bad='0.55', robust=True,
                    logscale=False, unit='min', show=True, **kwargs):

    # Define unit. Default is min
    if unit in ['second', 's', 'sec']:
        unitTime = 1
    else:
        unitTime = 60

    # Make dfUnqClusters indexed by annotation if not already
    if dfUnqClusters.index.name == 'annotation':
        dfUnqClusters2 = dfUnqClusters.copy()
    else:
        dfUnqClusters2 = dfUnqClusters.set_index('annotation')

    # See if user has provided an reference annotation; if so, compute
    # reference activity from reference annotation
    computeRefFromAnnt = isinstance(ref, basestring)

    # Define some commonly used transform functions
    if transform == 'kobs':
        def transformFunc(x): return 1. / x * unitTime
        logscale = False
        actLabel = r'$\mathrm{\mathsf{k_{obs}\ (min^{-1})}}$'
        if computeRefFromAnnt:
            ref = transformFunc(dfUnqClusters2.loc[ref][field])
        if ref:
            default_cmap = 'RdBu'
        else:
            default_cmap = 'YlOrRd'
    elif transform == 'logkobs':
        def transformFunc(x): return 1. / x * unitTime
        logscale = True
        actLabel = r'$\mathrm{\mathsf{k_{obs}\ (min^{-1})}}$'
        if computeRefFromAnnt:
            ref = transformFunc(dfUnqClusters2.loc[ref][field])
        if ref:
            default_cmap = 'RdBu'
        else:
            default_cmap = 'YlOrRd'
    elif transform == 'kobsfold':
        if computeRefFromAnnt:
            ref_kobs = (1. / dfUnqClusters2.loc[ref][field] * unitTime)
        def transformFunc(x): return ref_kobs / (1. / x * unitTime)
        logscale = False
        actLabel = r'$\mathrm{\mathsf{k_{obs}\ fold\ change}}$'
        ref = 1
        default_cmap = 'RdBu_r'
    elif transform == 'logkobsfold':
        if computeRefFromAnnt:
            ref_kobs = (1. / dfUnqClusters2.loc[ref][field] * unitTime)
        def transformFunc(x): return ref_kobs / (1. / x * unitTime)
        logscale = True
        actLabel = r'$\mathrm{\mathsf{k_{obs}\ fold\ change}}$'
        ref = 1
        default_cmap = 'RdBu_r'
    elif transform is None:
        def transformFunc(x): return x
        if computeRefFromAnnt:
            ref = transformFunc(dfUnqClusters2.loc[ref][field])
        if ref:
            default_cmap = 'RdBu'
        else:
            default_cmap = 'YlOrRd'
    else:
        transformFunc = transform
        if computeRefFromAnnt:
            ref = transformFunc(dfUnqClusters2.loc[ref][field])
        if ref:
            default_cmap = 'RdBu'
        else:
            default_cmap = 'YlOrRd'

    # Get the reordering index of 3 bases within the base triple if orderBaseTriple is provided
    # First base is the non-base-paired base, last two bases are base paired
    if orderBaseTriple is not None:
        listSeqPos = [mut[:-1] for mut in dfMut['mutations'][0]]
        orderListSeqPos = [listSeqPos.index(seqPos) for seqPos in orderBaseTriple]
    else:
        orderListSeqPos = [0, 1, 2]

    # Construct a new series of the transformed activity as specified by field
    # of the all the base triple variants
    seriesAct = transformFunc(dfUnqClusters2.loc[dfMut['annotation']][field])

    # Get multiindex from the list of tuples of mutations in dfMut
    # and reassign index as multiindex
    seriesAct.index = pd.MultiIndex.from_tuples(dfMut['mutations'])

    # Reorder index levels according to orderListSeqPos
    seriesAct = seriesAct.reorder_levels(orderListSeqPos)

    # Get the finite numbers for calculations
    seriesAct_finite = seriesAct[np.isfinite(seriesAct)]

    # Get all 4 mutations at the each of the base position
    listMutFirstBase = seriesAct.index.levels[0]
    listMutSecondBase = seriesAct.index.levels[1]
    listMutThirdBase = seriesAct.index.levels[2][::-1]

    # Find the multi-index of the WT variant
    WT_firstbase = [i for i, mut in enumerate(listMutFirstBase) if mut[0] == mut[-1]][0]
    WT_secondbase = [i for i, mut in enumerate(listMutSecondBase) if mut[0] == mut[-1]][0]
    WT_thirdbase = [i for i, mut in enumerate(listMutThirdBase) if mut[0] == mut[-1]][0]

    # Make figure
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    axes = axes.flatten()
    fig.tight_layout(rect=[0, 0, .88, 0.96])
    cbar_ax = fig.add_axes([.875, .1, .03, .8])

    # Set parameters for plotting
    # Set robust vmin and vmax
    if vmin is None:
        vmin = np.percentile(seriesAct_finite, 2) if robust else min(seriesAct_finite)
    if vmax is None:
        vmax = np.percentile(seriesAct_finite, 98) if robust else max(seriesAct_finite)
    # Take log if requested. Transform vmin and vmax if ref is provided
    if logscale:
        seriesAct = np.log10(seriesAct)
        vmin, vmax = np.log10(vmin), np.log10(vmax)
        if ref:
            ref = np.log10(ref)
            vlim = max(abs(vmin - ref), abs(vmax - ref))
            vmin, vmax = -vlim + ref, vlim + ref
            ticks = np.hstack([np.logspace(vmin, ref, 8), np.logspace(ref, vmax, 8)])
        else:
            ticks = np.logspace(vmin, vmax, 8)
        cbar_norm = mpl.colors.LogNorm(vmin=10**vmin, vmax=10**vmax)
        formatter = LogFormatter(10, labelOnlyBase=False)
    else:
        if ref:
            vlim = max(abs(vmin - ref), abs(vmax - ref))
            vmin, vmax = -vlim + ref, vlim + ref
        ticks = None
        cbar_norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        formatter = None
    # Define colormap
    if cmap is None:
        cmap = plt.get_cmap(default_cmap)
    elif isinstance(cmap, basestring):
        cmap = plt.get_cmap(cmap)
    cmap.set_bad(c_bad, 0.8)

    # Plotting the 4 different heatmaps
    for i, mut in enumerate(listMutFirstBase):

        # Get the 2D-ized DF of all the variants with the particular mutation in the first base
        currDF = seriesAct.loc[mut].unstack(level=-1).sort_index(axis=1, ascending=False)
        mask = ~np.isfinite(currDF)

        # Plot heatmap
        sns.heatmap(currDF, ax=axes[i], square=True, mask=mask, cbar=False,
                    vmin=vmin, vmax=vmax, cmap=cmap, center=ref)
        # Draw a box around the WT variant
        if i == WT_firstbase:
            _drawBox(axes[i], WT_thirdbase, 3-WT_secondbase)
        setproperties(ax=axes[i], yticklabelrot=90,
                      title=mut, titlefontsize=25,
                      suptitle=suptitle, suptitlefontsize=25,
                      borderwidth=0, tight=False, **kwargs)

        # Plot colorbar
        if i == 0:
            cbar = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap, norm=cbar_norm, ticks=ticks, format=formatter)
            cbar.ax.tick_params(labelsize=16)
            if actLabel is not None:
                cbar.set_label(actLabel, fontsize=22)

    if show:
        plt.show(block=False)

    return fig, axes, cbar_ax


# Convenient function to make base triple text
def notateBaseTriple(listSeqPos):
    return u"{} {} {}{}{}".format(listSeqPos[0], u'\u2022', listSeqPos[1], u'\u2013', listSeqPos[2])


# Rethink the complexity of wrapper functions, aka break plotBaseTriple into two functions
