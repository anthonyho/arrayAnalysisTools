# Anthony Ho, ahho@stanford.edu, 1/15/2015
# Last update 2/12/2016
"""Python module containing some handy plotting tools and common plot functions"""


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns
from scipy.stats import gaussian_kde
import seqlib



### ---------- Handy plotting tools ---------- ###


# Handy function to set the commonly used plot modifiers and apply to plot/figure
def setproperties(fig=None, ax=None, figsize=None,
                  suptitle=None, title=None,
                  legend=None, legendloc=1, legendwidth=2.5, legendbox=None,
                  xlabel=None, ylabel=None, xlim=None, ylim=None,
                  scix=False, sciy=False, scilimitsx=(-3, 3), scilimitsy=(-3, 3),
                  logx=False, logy=False, majorgrid=None, minorgrid=None,
                  borderwidth=2.5, tight=True, pad=1.6,
                  fontsize=None, legendfontsize=20, tickfontsize=20,
                  labelfontsize=20, titlefontsize=18, suptitlefontsize=20,
                  xticklabelrot=None, yticklabelrot=None, 
                  equal=False):
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
        legend = plt.legend(loc=legendloc, numpoints=1, fontsize=legendfontsize, frameon=legendbox)
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

    # Set equal aspect
    if equal:
        ax.set_aspect('equal', adjustable='box')


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



### ---------- Custom colormaps ---------- ###


# Define the RdYlBu_r2 colormap
# Essentially the same as RdYlBu_r but with the same darker shade
# of red and blue on the extremes as RdBu
def RdYlBu_r2():
    res = 10
    RdBu_r_cmap = cm.get_cmap("RdBu_r", res)
    RdBu_r_vals = RdBu_r_cmap(np.arange(res))
    
    RdYl_r_cmap = cm.get_cmap("RdYlBu_r", res)
    RdYl_r_vals = RdYl_r_cmap(np.arange(res))
    
    RdYl_r_vals[0] = RdBu_r_vals[0]
    RdYl_r_vals[-1] = RdBu_r_vals[-1]
    
    return mpl.colors.LinearSegmentedColormap.from_list("RdYlBu_r2", RdYl_r_vals)

# Handy function to show a colormap
def showColormap(cmap, numPoint=256):
    gradient = np.outer(np.ones(16), np.arange(numPoint))
    fig = plt.figure(figsize=(12, 2))
    ax = plt.imshow(gradient, cmap=cmap)
    ax.axes.get_yaxis().set_ticks([])



### ---------- Common plot functions ---------- ###


# Plot scatter plot colored by local density
def scatterDensity(data, data2=None, 
                   cmap=plt.cm.jet, colorbar=False, 
                   log=False, norm=None, sort=True, **kwargs):
    """Plot scatter plot colored by local density"""
    # Convert data to numpy array
    if data2 is None:
        data_np = np.array(data)
        x = data_np[:,0]
        y = data_np[:,1]
    else:
        x = np.array(data)
        y = np.array(data2)

    # Compute kernel density and sort data by z score if true
    xy = np.vstack([x,y])
    z = gaussian_kde(xy)(xy)
    if sort:
        idx = z.argsort()
        x, y, z = x[idx], y[idx], z[idx]

    # Plot and return
    if log and (norm is None):
        norm = LogNorm()
    ax1 = plt.scatter(x, y, c=z, cmap=cmap, norm=norm,
                      edgecolor='face', marker='o', **kwargs)
    if colorbar:
        ax2 = plt.colorbar()
        return ax1, ax2
    else:
        return ax1


# Generate the double mutant/cooperativity matrix
def doubleMutantMatrix(data, refVariant, libSeq, startPos, coop=False):
    """Auxiliary function to generate the doubel mutant matrix"""
    # Get library positions and create labels for mutants
    libPos = [i for (i, base) in enumerate(libSeq.upper()) if base == 'N']
    mutantLabels = [refVariant[i]+str(i+startPos)+otherBase 
                    for i in libPos for otherBase in seqlib.allOtherBases(refVariant[i])]
    
    # Grep the mutants and fill in the signals
    dim = len(mutantLabels)
    doubleMutantSignals = np.zeros(shape=(dim, dim))
    for i, mutation1 in enumerate(mutantLabels):
        for j, mutation2 in enumerate(mutantLabels):
            # Get the index for the degenerate base along the sequence
            pos1 = int(mutation1[1:-1]) - startPos
            otherBase1 = mutation1[-1]
            pos2 = int(mutation2[1:-1]) - startPos
            otherBase2 = mutation2[-1]
            # Create the current mutant sequence
            currSeq = list(refVariant)
            currSeq[pos1] = otherBase1
            currSeq[pos2] = otherBase2
            currSeq = ''.join(currSeq)
            # Grep the signal and fill in the double mutant matrix
            if not (currSeq in data.index):
                doubleMutantSignals[i, j] = np.nan
            elif pos1 != pos2:
                doubleMutantSignals[i, j] = data[currSeq]
            elif (pos1 == pos2) and (otherBase1 == otherBase2):
                doubleMutantSignals[i, j] = data[currSeq]
            else:
                doubleMutantSignals[i, j] = np.nan
    
    # Compute cooperativity if requested
    if coop:
        coopSignals = np.zeros(shape=(dim, dim))
        for i in xrange(dim):
            for j in xrange(dim):
                coopSignals[i, j] = doubleMutantSignals[i, i] + doubleMutantSignals[j, j] - doubleMutantSignals[i, j]
        doubleMutantSignals = coopSignals

    return doubleMutantSignals, mutantLabels


# Plot double mutant/cooperativity heatmap by grepping the activities of 
# all the double mutants given a reference sequence and a library sequence 
# with degenerate positions denoted with N
def doubleMutant(data, refVariant, libSeq, 
                 startPos=1, refSignal=None, normToRefSignal=True, coop=False,
                 vmin=None, vmax=None, cmap=None, center=0, cbarLabel=None,
                 triangle=None, invertY=True, linewidth=3, **kwargs):
    """Plot double mutant heatmap given a reference and library sequence"""
    # Define reference signal as the signal of the reference variant if 
    # refSignal not provided
    if refSignal is None:
        refSignal = data[refVariant]
    
    # Normalize data to reference signal if normToRefSignal=True
    if normToRefSignal:
        data_norm = data / refSignal
    else:
        data_norm = data

    # Generate the double mutant matrix
    doubleMutantSignals, mutantLabels = doubleMutantMatrix(data_norm, refVariant, 
                                                           libSeq, startPos, coop)

    # Create mask for triangular matrix if requested
    mask = np.zeros_like(doubleMutantSignals, dtype=bool)
    if triangle == 'lower':
        mask[np.tril_indices_from(mask)] = True
        mask = np.invert(mask)
    elif triangle == 'upper':
        mask[np.triu_indices_from(mask)] = True
        mask = np.invert(mask)

    # Plot the double mutant heatmap
    if cmap is None:
        cmap = RdYlBu_r2()
    ax = sns.heatmap(doubleMutantSignals, 
                     mask=mask, square=True, robust=True,
                     vmin=vmin, vmax=vmax, center=center, cmap=cmap, 
                     xticklabels=mutantLabels, yticklabels=mutantLabels, 
                     cbar_kws={'label': cbarLabel}, **kwargs)
    cax = plt.gcf().axes[-1]
    if invertY:
        ax.invert_yaxis()

    # Draw white lines separating the triplets
    dim = len(mutantLabels)
    for x in range(3, dim, 3):
        ax.plot([x, x], [0, dim], color='white', linewidth=linewidth)
    for y in range(3, dim, 3):
        ax.plot([0, dim], [y, y], color='white', linewidth=linewidth)

    return ax, cax
