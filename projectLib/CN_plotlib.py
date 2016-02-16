# Anthony Ho, ahho@stanford.edu, 1/26/2015
# Last update 2/15/2016
"""Python module containing plot functions for chemical nose project"""


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotlib


figureDir = 'figures/'


# Plot 
def plotSmQuenchScatter(data_aggPF, name, fullName):
    x = data_aggPF['QO_norm_signals']['median']
    y = data_aggPF['sm_norm_signals']['median']
    plt.figure(figsize=(7, 7))
    plotlib.scatterDensity(x, y, alpha=0.7, log=True)
    plotlib.setproperties(xlim=(-0.08, 1), ylim=(-0.08, 1), 
                          xlabel='Normalized quenched signals', 
                          ylabel='Normalized quenched signals \nin the prescence of '+name, 
                          title=fullName,
                          labelfontsize=25, tickfontsize=25, equal=True)
    plt.savefig(figureDir+'/'+name+'_smQuenchScatter.png')


# Plot
def plotSwitchQuenchScatter(data_aggPF, name, fullName):
    x = data_aggPF['QO_norm_signals']['median']
    y = data_aggPF['sm_norm_diff_signals_2']['median']
    plt.figure(figsize=(7, 7))
    plotlib.scatterDensity(x, y, alpha=0.7, log=True)
    plotlib.setproperties(xlim=(-0.08, 1), ylim=(-0.08, 1), 
                          xlabel='Normalized quenched signals', 
                          ylabel='Normalized switching signals', 
                          title=fullName,
                          labelfontsize=25, tickfontsize=25, equal=True)
    plt.savefig(figureDir+'/'+name+'_switchQuenchScatter.png')





# Plot raw signals histogram
def plotRawSignalsHist(data_aggPF, name, fullName):
    colors = sns.color_palette("Paired", 12)
    plt.figure(figsize=(10,7))
    ax = plt.subplot(1,1,1)
    sns.distplot(data_aggPF['sm_signals']['median'], color=colors[1], label='Cy3-read 7\' + QO1-BHQ1 + ' + fullName)
    sns.distplot(data_aggPF['QO_signals']['median'], color=colors[4], label='Cy3-read 7\' + QO1-BHQ1')
    sns.distplot(data_aggPF['r7_signals']['median'], color=colors[3], label='Cy3-read 7\' only')
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], loc=1, fontsize=16)  # Define legend and reverse order
    setproperties(xlim=(0, 800), 
                  xlabel='Median raw signal', ylabel='Relative frequency', 
                  labelfontsize=25, tickfontsize=25)
    plt.savefig('figures/'+name+'_rawSignalsHist.png')


# Plot normalized signals histogram
def plotNormSignalsHist(data_aggPF, name, fullName):
    colors = sns.color_palette("Paired", 12)
    plt.figure(figsize=(15,7))
    ax1 = plt.subplot(1,2,1)
    data_aggPF['sm_norm_signals']['median'].hist(bins=100, range=(0, 1),
                                                 color=colors[0], label='Cy3-read 7\' + QO1-BHQ1 + ' + fullName)
    data_aggPF['QO_norm_signals']['median'].hist(bins=100, range=(0, 1),
                                                 color=colors[1], label='Cy3-read 7\' + QO1-BHQ1')
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles[::-1], labels[::-1], loc=1, fontsize=16)  # Define legend and reverse order
    setproperties(xlabel='Median normalized signal', ylabel='Number of sequence variants', 
                  labelfontsize=25, tickfontsize=25)
    ax2 = plt.subplot(1,2,2)
    data_aggPF['sm_norm_signals']['median'].hist(bins=100, range=(0.1, 1),
                                                 color=colors[0], label='Cy3-read 7\' + QO1-BHQ1 + ' + fullName)
    data_aggPF['QO_norm_signals']['median'].hist(bins=100, range=(0.1, 1),
                                                 color=colors[1], label='Cy3-read 7\' + QO1-BHQ1')
    setproperties(xlabel='Median normalized signal', ylabel='Number of sequence variants', 
                  labelfontsize=25, tickfontsize=25)
    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend(handles[::-1], labels[::-1], loc=1, fontsize=16)  # Define legend and reverse order
    plt.savefig('figures/'+name+'_normSignalsHist.png')


# Plot the normalized signals along the sequence coordinate
def plotNormSignalsTrack(data_aggPF, name, fullName):
    colors = sns.color_palette("Paired", 12)
    plt.figure(figsize=(15,5))
    plt.plot(range(0, len(data_aggPF)), data_aggPF['sm_norm_signals']['median'], linewidth=0.5,
             color=colors[0], label='Cy3-read 7\' + QO1-BHQ1 + '+fullName)
    plt.plot(range(0, len(data_aggPF)), data_aggPF['QO_norm_signals']['median'], linewidth=0.5,
             color=colors[1], label='Cy3-read 7\' + QO1-BHQ1')
    setproperties(xlim=(0, len(data_aggPF)), ylim=(0, 1), 
                  xlabel='Sequence variants', ylabel='Meidan normalized signal', 
                  legend=True, legendloc=2, legendfontsize=16,
                  labelfontsize=25, tickfontsize=25)
    plt.savefig('figures/'+name+'_normSignalsTrack.png')


# Plot the normalized signals along the sequence coordinate, zoomed in to specific region
def plotNormSignalsTrackZoom(data_aggPF, name, fullName, start, end):
    colors = sns.color_palette("Paired", 12)
    plt.figure(figsize=(15,5))
    plt.plot(range(start, end), data_aggPF['sm_norm_signals']['median'][start:end], linewidth=1.2,
             color=colors[0], label='Cy3-read 7\' + QO1-BHQ1 + '+fullName)
    plt.plot(range(start, end), data_aggPF['QO_norm_signals']['median'][start:end], linewidth=1.2,
             color=colors[1], label='Cy3-read 7\' + QO1-BHQ1')
    setproperties(ylim=(0, 1),
                  xlabel='Sequence variants', ylabel='Meidan normalized signal', 
                  legend=True, legendloc=2, legendfontsize=16,
                  labelfontsize=25, tickfontsize=25)
    plt.savefig('figures/'+name+'_normSignalsTrackZoom_'+str(start)+'_'+str(end)+'.png')


# Plot the background substrated, normalized switching signals along the sequence coordinate
def plotNormSwitchingSignalsTrack(data_aggPF, name, fullName):
    colors = sns.color_palette("Paired", 12)
    plt.figure(figsize=(15,5))
    plt.plot(range(0, len(data_aggPF)), data_aggPF['sm_norm_diff_signals_2']['median'], linewidth=0.5,
             color=colors[5], label='Cy3-read 7\' + QO1-BHQ1 + ' + fullName + ', background substrated')
    setproperties(xlim=(0, len(data_aggPF)), ylim=(-0.2, 0.8), 
                  xlabel='Sequence variants', ylabel='Meidan normalized signal', 
                  legend=True, legendloc=2, legendfontsize=16,
                  labelfontsize=25, tickfontsize=25)
    plt.savefig('figures/'+name+'_normDiffSignalsTrack.png')


# Plot the background substrated, normalized switching signals along the sequence coordinate, zoomed in to specific region
def plotNormSwitchingSignalsTrackZoom(data_aggPF, name, fullName, start, end):
    colors = sns.color_palette("Paired", 12)
    plt.figure(figsize=(15,5))
    plt.plot(range(start, end), data_aggPF['sm_norm_diff_signals_2']['median'][start:end], linewidth=1.2,
             color=colors[5], label='Cy3-read 7\' + QO1-BHQ1 + ' + fullName + ', background substrated')
    setproperties(ylim=(-0.2, 0.8), 
                  xlabel='Sequence variants', ylabel='Meidan normalized signal', 
                  legend=True, legendloc=2, legendfontsize=16,
                  labelfontsize=25, tickfontsize=25)
    plt.savefig('figures/'+name+'_normDiffSignalsTrackZoom_'+str(start)+'_'+str(end)+'.png')


# Plot the background substrated, normalized switching signals along the sequence coordinate
def plotNormSwitchingSignalsSortedTrack(data_aggPF, name, fullName):
    colors = sns.color_palette("Paired", 12)
    plt.figure(figsize=(15,5))
    n = len(data_aggPF)
    x = range(0, n)
    sortedData = data_aggPF['sm_norm_diff_signals_2'].sort('median', ascending=False)
    y = sortedData['median']
    upper_bound = y + sortedData['sem']
    lower_bound = y - sortedData['sem']
    plt.plot(x, y, linewidth=2,
             color=colors[5], label='Cy3-read 7\' + QO1-BHQ1 + ' + fullName + ', background substrated, rank sorted')
    plt.fill_between(x, upper_bound, lower_bound, color=colors[4], alpha=0.5)
    setproperties(xlim=(0, n), ylim=(-0.2, 0.8), 
                  xlabel='Sorted sequence variants', ylabel='Meidan normalized signal', 
                  legend=True, legendloc=2, legendfontsize=16,
                  labelfontsize=25, tickfontsize=25)
    plt.savefig('figures/'+name+'_normDiffSignalsTrackSorted.png')


# Plot all normalized signals along the sequence coordinates
def plotNormSignalsSortedTrackAll(allData_aggPF, names, fullNames):
    colors = sns.color_palette("Paired", 12)
    plt.figure(figsize=(15,6))
    max_n = max([len(data_aggPF) for data_aggPF in allData_aggPF])
    for i, data_aggPF in enumerate(allData_aggPF):
        n = len(data_aggPF)
        x = range(0, n)
        plt.plot(x, data_aggPF['sm_norm_signals']['median'].order(ascending=False), linewidth=2, linestyle='dashed',
                 color=colors[i*2+1], label='Cy3-read 7\' + QO1-BHQ1 + '+fullNames[i])
        plt.plot(x, data_aggPF['QO_norm_signals']['median'].order(ascending=False), linewidth=2,
                 color=colors[i*2+1], label='Cy3-read 7\' + QO1-BHQ1')
    setproperties(xlim=(0, max_n), ylim=(0, 1), 
                  xlabel='Sorted sequence variants', ylabel='Meidan normalized signal', 
                  legend=True, legendloc=1, legendfontsize=16,
                  labelfontsize=25, tickfontsize=25)
    plt.savefig('figures/all6_normSignalsTrackSorted.png')


# Plot all background substrated, normalized switching signals along the sequence coordinates
def plotNormSwitchingSignalsSortedTrackAll(allData_aggPF, names, fullNames):
    colors = sns.color_palette("Paired", 12)
    plt.figure(figsize=(15,5))
    max_n = max([len(data_aggPF) for data_aggPF in allData_aggPF])
    for i, data_aggPF in enumerate(allData_aggPF):
        n = len(data_aggPF)
        x = range(0, n)
        sortedData = data_aggPF['sm_norm_diff_signals_2'].sort('median', ascending=False)
        y = sortedData['median']
        #upper_bound = y + sortedData['sem']
        #lower_bound = y - sortedData['sem']
        plt.plot(x, y, linewidth=2,
                 color=colors[i*2+1], label='Cy3-read 7\' + QO1-BHQ1 + ' + fullNames[i] + ', background substrated, rank sorted')
        #plt.fill_between(x, upper_bound, lower_bound, color=colors[i*2], alpha=0.5)
    setproperties(xlim=(0, max_n), ylim=(-0.2, 0.8), 
                  xlabel='Sorted sequence variants', ylabel='Meidan normalized signal', 
                  legend=True, legendloc=1, legendfontsize=16,
                  labelfontsize=25, tickfontsize=25)
    plt.savefig('figures/all6_normDiffSignalsTrackSorted.png')

