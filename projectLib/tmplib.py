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



# Plot box and/or violin plots of variant activities along with their counts
# This is the most low level of the plotActCount-related script that is not
# intended to be used on its own; consider using plotVariants or
# plotSingleMutants instead
def plotActCount(listVarDF, field, listName=None, color=None,
                 transform=None, ref=None, plotmode='b', bootstrap=2000, inner=None,
                 logAct=False, logCount=True,
                 actLabel=None, actLim=None,
                 barwidth=0.75, xticklabelrot=None,
                 figsize=None, figunitheight=4, figunitwidth=2, maxfigwidth=32,
                 _xticks=None, show=True, detLim=40, detLimRef=1, **kwargs):
    """Plot box and/or violin plots of variant activities along with their counts

    Required arguments:
      listVarDF - a list of dataframe, each of them belongs to a single variant
      field - the name of the column of the dataframes to be plotted
    """
    # Define parameters
    numVar = len(listVarDF)
    numSubplots = 1  # np.sum(['b' in plotmode, 'v' in plotmode]) + 1
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
    if actLim is not None:
        logActLim = np.log10(actLim)
    else:
        logActLim = actLim

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
        figheight = 6 * numSubplots
        figsize = (figwidth, figheight)
#    fig, axes = plt.subplots(numSubplots, 1, sharex=False, figsize=figsize)
    fig = plt.figure(figsize=figsize)
    axBox = plt.gca()
    fig.patch.set_facecolor('w')
#    axCount = axes[-1]
#    if 'b' in plotmode and 'v' in plotmode:
#        axBox = axes[0]
#        axViolin = axes[1]
#    elif 'b' in plotmode and 'v' not in plotmode:
#        axBox = axes[0]
#    elif 'b' not in plotmode and 'v' in plotmode:
#        axViolin = axes[0]    

    # Make box plot
    if 'b' in plotmode:
        if ref:
            axBox.plot(xlim, [ref]*2, linestyle='--', color='k', linewidth=2)
            axBox.plot(xlim, [detLim/detLimRef]*2, linestyle='--', color='m', linewidth=2)
        listAct2 = np.array([np.median(act) for act in listAct])
        listErr = np.array([1.253*np.std(act)/np.sqrt(len(act)) for act in listAct])
        axBox.bar(_xticks-barwidth/2, listAct2, color=color, width=barwidth,
                  yerr=listErr, capsize=5, ecolor='k', error_kw={'capthick': 1.5, 'elinewidth': 1.5})
#        sns.boxplot(listAct, positions=_xticks, color=color,
#                    showmeans=True, bootstrap=bootstrap,
#                    widths=barwidth, ax=axBox)
        axBox.set_xticks(_xticks)
        axBox.set_xticklabels(listName)
        setproperties(ax=axBox, logy=logAct, ylabel=actLabel, ylim=actLim,
                      majorgrid=True, minorgrid=logAct,
                      xticklabelrot=xticklabelrot, xlim=xlim, pad=2.5, **kwargs)

    # Make count bar chart
#    axCount.bar(_xticks-barwidth/2, listFilteredCount, width=barwidth)
#    axCount.bar(_xticks-barwidth/2, listRawCount-listFilteredCount, width=barwidth,
#                bottom=listFilteredCount, color=sns.xkcd_rgb["tangerine"])
#    axCount.set_xticklabels(listName)
#    setproperties(ax=axCount, logy=logCount, ylabel='Number of clusters',
#                  ylim=(1, None), xlim=xlim, xticklabelrot=xticklabelrot, **kwargs)

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
            wtref = (1. / df2.loc[['WT:0:0:0::::']][field] * unitTime).median()
        def transformFunc(x): return ref / (1. / x * unitTime)
        kwargs['logAct'] = True
        kwargs['actLabel'] = r'$\mathrm{\mathsf{k_{obs}\ fold\ change}}$'
        kwargs['ref'] = 1
        kwargs['detLimRef'] = wtref/ref
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
            mNames = [name[:-1] if i % 3 == 1 else ' ' for i, name in enumerate(mNames)]
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
