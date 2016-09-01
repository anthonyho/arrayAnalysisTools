# Anthony Ho, ahho@stanford.edu, 6/4/2016
# Last update 6/4/2016
"""Python module containing plot functions for chemical nose project"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import plotlib


figureDir = 'figures/'




### ---------- Scatter plots of raw and normalized signals ---------- ###


# Plot normalized switching signals against normalized quenched signals
def plotScatterNorm_QO_sm(data_aggPF, name, fullName):
    x = data_aggPF['quenched_norm']['median']
    y = data_aggPF['switch2_norm']['median']
    plt.figure(figsize=(7, 7))
    plotlib.scatterDensity(x, y, alpha=0.7, log=True)
    plotlib.setproperties(xlim=(-0.09, 1), ylim=(-0.09, 1), 
                          xlabel='Normalized quenched signals', 
                          ylabel='Normalized switching signals', 
                          title=fullName,
                          labelfontsize=25, tickfontsize=25, equal=True)
    plt.savefig(figureDir+'/'+name+'_scatterNorm_QO_sm.png')


# Plot normalized switching signals against normalized quenched signals
def plotScatterNorm_QO_sm_linear(data_aggPF, name, fullName):
    x = data_aggPF['quenched_norm']['median']
    y = data_aggPF['switch2_norm']['median']
    plt.figure(figsize=(7, 7))
    plotlib.scatterDensity(x, y, alpha=0.7, log=False)
    plotlib.setproperties(xlim=(-0.09, 1), ylim=(-0.09, 1), 
                          xlabel='Normalized quenched signals', 
                          ylabel='Normalized switching signals', 
                          title=fullName,
                          labelfontsize=25, tickfontsize=25, equal=True)
    plt.savefig(figureDir+'/'+name+'_scatterNorm_QO_sm_linear.png')
