# Anthony Ho, ahho@stanford.edu, 1/10/2016
# Last update 1/10/2015
"""Library containing plotting functions"""

import matplotlib.pyplot as plt
import seaborn as sns
import plotlib
import CN_globalVars


def plotPredictedConcMatrixPureSamples(predictedConcMatrix, setup, nAptamers=None, 
                                       cat_colors=None, vmax=350, figsize=(9.15, 9), fig_dir=None):

    # Define cat_colors
    if cat_colors is None:
        listCatColors = [CN_globalVars.catcolors[sm] for sm in predictedConcMatrix.columns]
    elif cat_colors is False:
        listCatColors = None
    else:
        listCatColors = cat_colors

    # Make plot
    cg = sns.clustermap(predictedConcMatrix,
                        row_colors=listCatColors, col_colors=listCatColors,
                        figsize=figsize, row_cluster=False, col_cluster=False, vmax=vmax, vmin=0,
                        cbar_kws={'label': 'Predicted concentration (uM)'})
    cax = plt.gcf().axes[-1]
    cax.set_position([0.13, .2, .03, .45])

    # Set miscel properties
    if cat_colors is False:
        title_y = 1.03
    else: 
        title_y = 1.1
    if nAptamers is None:
        cg.ax_heatmap.set_title(setup, y=title_y)
    else:
        cg.ax_heatmap.set_title(setup+', n='+str(nAptamers), y=title_y)
    plotlib.setproperties(ax=cg.ax_heatmap, tight=False, 
                          xlabel='Ligand standards', ylabel='Ligand predictions',
                          fontsize=18, xticklabelrot=90, yticklabelrot=0)
    plotlib.setproperties(ax=cax, tight=False, fontsize=18, yticklabelrot=0)

    # Save figure
    if fig_dir is not None:
        cg.savefig(fig_dir+'/pureSample_conMat_'+setup+'.png')
        cg.savefig(fig_dir+'/pureSample_conMat_'+setup+'.eps')

