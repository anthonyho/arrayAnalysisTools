# Anthony Ho, ahho@stanford.edu, 1/10/2016
# Last update 1/10/2015
"""Library containing plotting functions"""

import matplotlib.pyplot as plt
import seaborn as sns
import plotlib


def plotPredictedConcMatrix(predictedConcMatrix, setup, nAptamers=None, cat_colors=None, 
                            vmax=350, figsize=(8, 7), fig_dir=None):
    plt.figure(figsize=(8, 7))
    cg = sns.clustermap(predictedConcMatrix,
                        row_colors=cat_colors, col_colors=cat_colors,
                        row_cluster=False, col_cluster=False, vmax=vmax, vmin=0,
                        cbar_kws={'label': 'Predicted concentration (uM)'})
    cax = plt.gcf().axes[-1]
    cax.set_position([0.13, .2, .03, .45])
    if cat_colors is None:
        title_y = 1.05
    else: 
        title_y = 1.1
    if nAptamers is None:
        cg.ax_heatmap.set_title(setup, y=title_y)
    else:
        cg.ax_heatmap.set_title(setup+', n='+str(nAptamers), y=title_y)
    plotlib.setproperties(ax=cg.ax_heatmap, tight=False, 
                          xlabel='Pure sample of ...', ylabel='Predicted as ...',
                          fontsize=18, xticklabelrot=90, yticklabelrot=0)
    plotlib.setproperties(ax=cax, tight=False, fontsize=18, yticklabelrot=0)
    if fig_dir is not None:
        cg.savefig(fig_dir+'/pureSample_conMat_'+setup+'.png')
        cg.savefig(fig_dir+'/pureSample_conMat_'+setup+'.eps')
