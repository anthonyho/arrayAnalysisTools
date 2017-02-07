# Anthony Ho, ahho@stanford.edu, 1/10/2016
# Last update 2/6/2015
"""Library containing plotting functions"""


import matplotlib.pyplot as plt
import seaborn as sns
import plotlib
import CN_globalVars


def plotPredictedConcMatrix(fitResults, setup,
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
        listCatColorsRow = [CN_globalVars.catcolors[sample] for sample in matrix.index]
    elif catColorsRow is False:
        listCatColorsRow = None
    else:
        listCatColorsRow = catColorsRow

    if catColorsCol is None:
        try:
            listCatColorsCol = [CN_globalVars.catcolors[sample] for sample in matrix.columns]
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
