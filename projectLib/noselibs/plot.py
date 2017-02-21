# Anthony Ho, ahho@stanford.edu, 1/10/2016
# Last update 2/7/2015
"""Library containing plotting functions"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import seaborn as sns
import plotlib
import globalvars


colors = sns.color_palette("Paired", 12)


def plotPredictedConcMatrix(fitResults, setup='',
                            nAptamers=None, vmax=None,
                            catColorsRow=None, catColorsCol=None, 
                            figsize=(9.1, 9), fig_dir=None):

    # Unpack object
    matrix = fitResults.predictedConcMatrix
    if nAptamers is None:
        nAptamers = fitResults.results[matrix.columns[0]].ndata
    try:
        vmax = 1.5 * fitResults.currConc
    except TypeError:
        pass

    # Define catColors
    if catColorsRow is None:
        listCatColorsRow = [globalvars.ligand_catcolors[ligand] for ligand in matrix.index]
    elif catColorsRow is False:
        listCatColorsRow = None
    else:
        listCatColorsRow = catColorsRow

    if catColorsCol is None:
        try:
            listCatColorsCol = [globalvars.ligand_catcolors[sample] for sample in matrix.columns]
        except KeyError:
            pass
    elif catColorsCol is False:
        listCatColorsCol = None
    else:
        listCatColorsCol = catColorsCol

    # Make plot
    cg = sns.clustermap(matrix,
                        row_colors=listCatColorsRow, col_colors=listCatColorsCol,
                        figsize=figsize, row_cluster=False, col_cluster=False, vmax=vmax, vmin=0,
                        cbar_kws={'label': 'Predicted concentration (uM)'})
    cax = plt.gcf().axes[-1]
    cax.set_position([0.13, .2, .03, .45])

    # Set miscel properties
    if catColorsCol is False:
        title_y = 1.03
    else: 
        title_y = 1.1
    if nAptamers:
        cg.ax_heatmap.set_title(setup+', n='+str(nAptamers), y=title_y)
    else:
        cg.ax_heatmap.set_title(setup, y=title_y)
    plotlib.setproperties(ax=cg.ax_heatmap, tight=False, 
                          xlabel='Samples', ylabel='Predictions',
                          fontsize=18, xticklabelrot=90, yticklabelrot=0)
    plotlib.setproperties(ax=cax, tight=False, fontsize=18, yticklabelrot=0)

    # Save figure
    if fig_dir is not None:
        cg.savefig(fig_dir+'/predictConMat_'+setup+'.png')
        cg.savefig(fig_dir+'/predictConMat_'+setup+'.eps')

    return cg


def plotFitStatus(fitResults, setup='',
                  metric='IERMSLE', figsize=(14, 8), fig_dir=None):

    fig = plt.figure(figsize=figsize)

    # Compute log2 fold change of fitted vs true unweighted chi2
    fittedRedChi = pd.Series({sample: fitResults.reportFit(sample, 'redchi', weighted=False, params='fitted') 
                            for sample in fitResults.listSamples})
    trueRedChi = pd.Series({sample: fitResults.reportFit(sample, 'redchi', weighted=False, params='true') 
                            for sample in fitResults.listSamples})
    redChiFoldChange = np.log2(fittedRedChi / trueRedChi)
    # Plot 
    ax1 = plt.subplot(3, 2, 1)
    sns.barplot(x=fitResults.listSamples, y=redChiFoldChange, 
                color=colors[3], edgecolor=colors[3])
    plotlib.setproperties(ax=ax1, fontsize=16, tight=False,
                          ylabel='Log2 fold change\nof fitted vs true\nunweighted chi2')
    plt.setp(ax1.get_xticklabels(), visible=False)

    # Compute log2 fold change of fitted vs true weighted chi2
    fittedRedChi = pd.Series({sample: fitResults.reportFit(sample, 'redchi', weighted=True, params='fitted') 
                            for sample in fitResults.listSamples})
    trueRedChi = pd.Series({sample: fitResults.reportFit(sample, 'redchi', weighted=True, params='true')
                            for sample in fitResults.listSamples})
    redChiFoldChange = np.log2(fittedRedChi / trueRedChi)
    # Plot
    ax2 = plt.subplot(3, 2, 3, sharex=ax1)
    sns.barplot(x=fitResults.listSamples, y=redChiFoldChange, 
                color=colors[3], edgecolor=colors[3])
    plotlib.setproperties(ax=ax2, fontsize=16, tight=False,
                          ylabel='Log2 fold change\nof fitted vs true\nweighted chi2')
    plt.setp(ax2.get_xticklabels(), visible=False)

    # Compute performance metric
    performance = fitResults.evaluatePerformance(metric)

    # Plot performance metric
    ax3 = plt.subplot(3, 2, 5, sharex=ax1)
    sns.barplot(x=fitResults.listSamples, y=performance, 
                color=colors[1], edgecolor=colors[3])
    plotlib.setproperties(ax=ax3, xlabel='Samples', ylabel=metric,
                          tight=False,
                          xticklabelrot=90, fontsize=16)

    # Compute residual matrics
    fittedRedMat = pd.DataFrame({sample: fitResults.reportFit(sample, 'residual', weighted=True, params='fitted') 
                                 for sample in fitResults.listSamples})
    trueRedMat = pd.DataFrame({sample: fitResults.reportFit(sample, 'residual', weighted=True, params='true') 
                               for sample in fitResults.listSamples})

    # Plot distribution of mean residuals across ligands of all aptamers
    ax4 = plt.subplot(3, 2, 2)
    sns.distplot(fittedRedMat.mean(axis=1), label='fitted')
    sns.distplot(trueRedMat.mean(axis=1), label='true')
    xlim = ax4.get_xlim()
    x = np.linspace(xlim[0], xlim[1], 1000)
    plt.plot(x, mlab.normpdf(x, 0, 1./np.sqrt(len(fitResults.results))), label='theoretical')
    plotlib.setproperties(ax=ax4, fontsize=16, tight=False, 
                          xlabel='Weighted residuals', ylabel='Density', 
                          legend=True, legendfontsize=12)

    plt.suptitle(setup, y = 1.01, fontsize=16)
    fig.tight_layout()

    # Save figure
    if fig_dir is not None:
        plt.savefig(fig_dir+'/fitStatus_'+setup+'.png')
        plt.savefig(fig_dir+'/fitStatus_'+setup+'.eps')
