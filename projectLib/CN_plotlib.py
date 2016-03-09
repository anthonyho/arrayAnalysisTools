# Anthony Ho, ahho@stanford.edu, 1/26/2015
# Last update 3/4/2016
"""Python module containing plot functions for chemical nose project"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import plotlib


figureDir = 'figures/'




### ---------- Scatter plots of raw and normalized signals ---------- ###


# Plot raw switching signals against raw all-cluster signals
def plotScatter_r7_sm(data_aggPF, name, fullName):
    x = data_aggPF['r7_signals']['median']
    y = data_aggPF['sm_signals']['median']
    plt.figure(figsize=(7, 7))
    plotlib.scatterDensity(x, y, alpha=0.7, log=True)
    ub = max(np.percentile(x, 99.9), np.percentile(y, 99.9)) * 1.1
    plotlib.setproperties(xlim=(0, ub), ylim=(0, ub), 
                          xlabel='All cluster signals', 
                          ylabel='Switching signals', 
                          title=fullName,
                          labelfontsize=25, tickfontsize=25, equal=True)
    plt.savefig(figureDir+'/'+name+'_scatter_r7_sm.png')


# Plot raw quenched signals against raw all-cluster signals 
def plotScatter_r7_QO(data_aggPF, name, fullName):
    x = data_aggPF['r7_signals']['median']
    y = data_aggPF['QO_signals']['median']
    plt.figure(figsize=(7, 7))
    plotlib.scatterDensity(x, y, alpha=0.7, log=True)
    ub = max(np.percentile(x, 99.9), np.percentile(y, 99.9)) * 1.1
    plotlib.setproperties(xlim=(0, ub), ylim=(0, ub), 
                          xlabel='All cluster signals', 
                          ylabel='Quenched signals', 
                          title=fullName,
                          labelfontsize=25, tickfontsize=25, equal=True)
    plt.savefig(figureDir+'/'+name+'_scatter_r7_QO.png')


# Plot raw switching signals against raw quenched signals
def plotScatter_QO_sm(data_aggPF, name, fullName):
    x = data_aggPF['QO_signals']['median']
    y = data_aggPF['sm_signals']['median']
    plt.figure(figsize=(7, 7))
    plotlib.scatterDensity(x, y, alpha=0.7, log=True)
    ub = max(np.percentile(x, 99.9), np.percentile(y, 99.9)) * 1.1
    plotlib.setproperties(xlim=(0, ub), ylim=(0, ub), 
                          xlabel='Quenched signals', 
                          ylabel='Switching signals', 
                          title=fullName,
                          labelfontsize=25, tickfontsize=25, equal=True)
    plt.savefig(figureDir+'/'+name+'_scatter_QO_sm.png')


# Plot normalized switching signals against normalized quenched signals
def plotScatterNorm_QO_sm(data_aggPF, name, fullName):
    x = data_aggPF['QO_norm_signals']['median']
    y = data_aggPF['sm_norm_signals']['median']
    plt.figure(figsize=(7, 7))
    plotlib.scatterDensity(x, y, alpha=0.7, log=True)
    plotlib.setproperties(xlim=(-0.09, 1), ylim=(-0.09, 1), 
                          xlabel='Normalized quenched signals', 
                          ylabel='Normalized switching signals', 
                          title=fullName,
                          labelfontsize=25, tickfontsize=25, equal=True)
    plt.savefig(figureDir+'/'+name+'_scatterNorm_QO_sm.png')


# Plot normalized differential signals 1 against normalized quenched signals
def plotScatterNorm_QO_ds1(data_aggPF, name, fullName):
    x = data_aggPF['QO_norm_signals']['median']
    y = data_aggPF['sm_norm_diff_signals_1']['median']
    plt.figure(figsize=(7, 7))
    plotlib.scatterDensity(x, y, alpha=0.7, log=True)
    plotlib.setproperties(xlim=(-0.09, 1), ylim=(-0.09, 1), 
                          xlabel='Normalized quenched signals', 
                          ylabel='Normalized differential signals 1', 
                          title=fullName,
                          labelfontsize=25, tickfontsize=25, equal=True)
    plt.savefig(figureDir+'/'+name+'_scatterNorm_QO_ds1.png')


# Plot normalized differential signals 2 against normalized quenched signals
def plotScatterNorm_QO_ds2(data_aggPF, name, fullName):
    x = data_aggPF['QO_norm_signals']['median']
    y = data_aggPF['sm_norm_diff_signals_2']['median']
    plt.figure(figsize=(7, 7))
    plotlib.scatterDensity(x, y, alpha=0.7, log=True)
    plotlib.setproperties(xlim=(-0.09, 1), ylim=(-0.09, 1), 
                          xlabel='Normalized quenched signals', 
                          ylabel='Normalized differential signals 2', 
                          title=fullName,
                          labelfontsize=25, tickfontsize=25, equal=True)
    plt.savefig(figureDir+'/'+name+'_scatterNorm_QO_ds2.png')


# Plot normalized differential signals 1 against normalized switching signals
def plotScatterNorm_sm_ds1(data_aggPF, name, fullName):
    x = data_aggPF['sm_norm_signals']['median']
    y = data_aggPF['sm_norm_diff_signals_1']['median']
    plt.figure(figsize=(7, 7))
    plotlib.scatterDensity(x, y, alpha=0.7, log=True)
    plotlib.setproperties(xlim=(-0.09, 1), ylim=(-0.09, 1), 
                          xlabel='Normalized switching signals', 
                          ylabel='Normalized differential signals 1', 
                          title=fullName,
                          labelfontsize=25, tickfontsize=25, equal=True)
    plt.savefig(figureDir+'/'+name+'_scatterNorm_sm_ds1.png')


# Plot normalized differential signals 2 against normalized switching signals
def plotScatterNorm_sm_ds2(data_aggPF, name, fullName):
    x = data_aggPF['sm_norm_signals']['median']
    y = data_aggPF['sm_norm_diff_signals_2']['median']
    plt.figure(figsize=(7, 7))
    plotlib.scatterDensity(x, y, alpha=0.7, log=True)
    plotlib.setproperties(xlim=(-0.09, 1), ylim=(-0.09, 1), 
                          xlabel='Normalized switching signals', 
                          ylabel='Normalized differential signals 2', 
                          title=fullName,
                          labelfontsize=25, tickfontsize=25, equal=True)
    plt.savefig(figureDir+'/'+name+'_scatterNorm_sm_ds2.png')


# Plot normalized differential signals 2 against normalized differential signals 1
def plotScatterNorm_ds1_ds2(data_aggPF, name, fullName):
    x = data_aggPF['sm_norm_diff_signals_1']['median']
    y = data_aggPF['sm_norm_diff_signals_2']['median']
    plt.figure(figsize=(7, 7))
    plotlib.scatterDensity(x, y, alpha=0.7, log=True)
    plotlib.setproperties(xlim=(-0.09, 1), ylim=(-0.09, 1), 
                          xlabel='Normalized differential signals 1', 
                          ylabel='Normalized differential signals 2', 
                          title=fullName,
                          labelfontsize=25, tickfontsize=25, equal=True)
    plt.savefig(figureDir+'/'+name+'_scatterNorm_ds1_ds2.png')


# Plot bootstrapped error vs median for normalized differential signals 2
def plotScatterNorm_ds2_bserr(data_aggPF, name, fullName):
    x = data_aggPF['sm_norm_diff_signals_2']['median']
    y = data_aggPF['sm_norm_diff_signals_2']['bserr'].apply(lambda x: float(x.split(',')[1][:-1]) - float(x.split(',')[0][1:]))
    plt.figure(figsize=(7, 7))
    plotlib.scatterDensity(x, y, alpha=0.7, log=True)
    plotlib.setproperties(xlabel='Normalized differential signals', 
                          ylabel='Bootstrapped error', 
                          title=fullName,
                          labelfontsize=25, tickfontsize=25, equal=True)
    plt.savefig(figureDir+'/'+name+'_scatterNorm_ds2_bserr.png')


# Plot fractional bootstrapped error vs median for normalized differential signals 2
def plotScatterNorm_ds2_fbserr(data_aggPF, name, fullName):
    x = data_aggPF['sm_norm_diff_signals_2']['median']
    y = data_aggPF['sm_norm_diff_signals_2']['bserr'].apply(lambda x: float(x.split(',')[1][:-1]) - float(x.split(',')[0][1:])) 
    df = pd.concat([x, y / x], axis=1).replace([np.inf, -np.inf], np.nan).dropna()
    plt.figure(figsize=(7, 7))
    plotlib.scatterDensity(df[0], df[1], alpha=0.7, log=True)
    plotlib.setproperties(xlabel='Normalized differential signals',
                          ylabel='Fractional bootstrapped error', 
                          title=fullName,
                          labelfontsize=25, tickfontsize=25, equal=True)
    plt.savefig(figureDir+'/'+name+'_scatterNorm_ds2_fbserr.png')




### ---------- Histograms of raw and normalized signals ---------- ###


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




### ---------- Signal tracks  ---------- ###


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




### ---------- Double mutant/cooperativity heatmaps  ---------- ###


# Plot double mutant/cooperativity heatmap for one small molecule
def doubleMutant0(allData_aggPF, mutantSM, dataSM, names, coop=False):
    # Define variables
    data = allData_aggPF[dataSM]['sm_norm_diff_signals_2']['median']
    mutantRefVariant = allData_aggPF[mutantSM]['sm_norm_diff_signals_2']['median'].idxmax()
    dataRefVariant = allData_aggPF[dataSM]['sm_norm_diff_signals_2']['median'].idxmax()
    refSignal = allData_aggPF[dataSM]['sm_norm_diff_signals_2']['median'][dataRefVariant]
    libSeq = 'NNGGATTTTCCNNNNACGAAGTNNTCCCGAG'
    startPos = 14
    cbarLabel = 'Switching signals normalized to '+names[dataSM].lower()+'0'
    title = 'Mutants of '+names[mutantSM].lower()+'0, in the presence of '+names[dataSM]
    
    # Compute vmin and vmax manually (assuming robust=True) so that the color bar
    # scale is the same across the same dataRefVariant
    normRefData = allData_aggPF[dataSM]['sm_norm_diff_signals_2']['median'] / refSignal
    doubleMutantSignals, mutantLabels = plotlib.doubleMutantMatrix(normRefData, dataRefVariant, libSeq, startPos, coop)
    doubleMutantSignals = doubleMutantSignals[~np.isnan(doubleMutantSignals)]
    vmin = np.percentile(doubleMutantSignals, 1)
    vmax = np.percentile(doubleMutantSignals, 99)
    vlim = max(abs(vmin), abs(vmax))
    vmin, vmax = -vlim, vlim
    
    # Make heatmap
    plt.figure(figsize=(12,10))
    ax, cax = plotlib.doubleMutant(data, mutantRefVariant, libSeq, 
                                   startPos=startPos, refSignal=refSignal, coop=coop,
                                   vmin=vmin, vmax=vmax, cbarLabel=cbarLabel,
                                   triangle='lower', invertY=False, linewidth=3)
    cax.tick_params(labelsize=20)
    cax.yaxis.label.set_fontsize(20)
    plotlib.setproperties(title=title, labelfontsize=20, titlefontsize=20)
    
    # Save to file
    if coop:
        plt.savefig(figureDir+'/'+names[dataSM]+'_'+names[mutantSM].lower()+'0_coop.png')
    else:
        plt.savefig(figureDir+'/'+names[dataSM]+'_'+names[mutantSM].lower()+'0_doubleMutant.png')
        
    return ax, cax


# Plot double mutant/cooperativity heatmap for one small molecule
def doubleMutant1(allData_aggPF, mutantSM, targetSM, normSM, names, coop=False):
    # Define variables
    data = allData_aggPF[targetSM]['sm_norm_diff_signals_2']['median']
    mutantRefVariant = allData_aggPF[mutantSM]['sm_norm_diff_signals_2']['median'].idxmax()
    targetRefVariant = allData_aggPF[targetSM]['sm_norm_diff_signals_2']['median'].idxmax()
    normRefVariant = allData_aggPF[normSM]['sm_norm_diff_signals_2']['median'].idxmax()
    targetRefSignal = allData_aggPF[targetSM]['sm_norm_diff_signals_2']['median'][targetRefVariant]
    normRefSignal = allData_aggPF[normSM]['sm_norm_diff_signals_2']['median'][normRefVariant]
    libSeq = 'NNGGATTTTCCNNNNACGAAGTNNTCCCGAG'
    startPos = 14
    cbarLabel = 'Switching signals normalized to '+names[targetSM].lower()+'0'
    title = 'Mutants of '+names[mutantSM].lower()+'0, in the presence of '+names[targetSM]
    
    # Compute vmin and vmax manually (assuming robust=True) so that the color bar
    # scale is the same across the same dataRefVariant
    normRefData = allData_aggPF[normSM]['sm_norm_diff_signals_2']['median'] / normRefSignal
    doubleMutantSignals, mutantLabels = plotlib.doubleMutantMatrix(normRefData, normRefVariant, libSeq, startPos, coop)
    doubleMutantSignals = doubleMutantSignals[~np.isnan(doubleMutantSignals)]
    vmin = np.percentile(doubleMutantSignals, 1)
    vmax = np.percentile(doubleMutantSignals, 99)
    vlim = max(abs(vmin), abs(vmax))
    vmin, vmax = -vlim, vlim
    
    # Make heatmap
    plt.figure(figsize=(12,10))
    ax, cax = plotlib.doubleMutant(data, mutantRefVariant, libSeq, 
                                   startPos=startPos, refSignal=targetRefSignal, coop=coop,
                                   vmin=vmin, vmax=vmax, cbarLabel=cbarLabel,
                                   triangle='lower', invertY=False, linewidth=3)
    cax.tick_params(labelsize=20)
    cax.yaxis.label.set_fontsize(20)
    plotlib.setproperties(title=title, labelfontsize=20, titlefontsize=20)
    
    # Save to file
    if coop:
        plt.savefig(figureDir+'/'+names[mutantSM].lower()+'0_'+names[targetSM]+'_'+names[normSM].lower()+'0_coop.png')
    else:
        plt.savefig(figureDir+'/'+names[mutantSM].lower()+'0_'+names[targetSM]+'_'+names[normSM].lower()+'0_doubleMutant.png')
        
    return ax, cax


# Plot double mutant/cooperativity heatmap for one small molecule
def doubleMutant(allData_aggPF, mutantSM, targetSM, normSM, names, norm=False, coop=False):
    # Define variables
    data = allData_aggPF[targetSM]['sm_norm_diff_signals_2']['median']
    mutantRefVariant = allData_aggPF[mutantSM]['sm_norm_diff_signals_2']['median'].idxmax()
    targetRefVariant = allData_aggPF[targetSM]['sm_norm_diff_signals_2']['median'].idxmax()
    normRefVariant = allData_aggPF[normSM]['sm_norm_diff_signals_2']['median'].idxmax()
    if norm:
        targetRefSignal = allData_aggPF[targetSM]['sm_norm_diff_signals_2']['median'][targetRefVariant]
        normRefSignal = allData_aggPF[normSM]['sm_norm_diff_signals_2']['median'][normRefVariant]
    else:
        targetRefSignal = 1
        normRefSignal = 1
    libSeq = 'NNGGATTTTCCNNNNACGAAGTNNTCCCGAG'
    startPos = 14
    cbarLabel = 'Switching signals normalized to '+names[targetSM].lower()+'0'
    title = 'Mutants of '+names[mutantSM].lower()+'0, in the presence of '+names[targetSM]
    
    # Compute vmin and vmax manually (assuming robust=True) so that the color bar
    # scale is the same across the same dataRefVariant
    normRefData = allData_aggPF[normSM]['sm_norm_diff_signals_2']['median'] / normRefSignal
    doubleMutantSignals, mutantLabels = plotlib.doubleMutantMatrix(normRefData, normRefVariant, libSeq, startPos, coop)
    doubleMutantSignals = doubleMutantSignals[~np.isnan(doubleMutantSignals)]
    vmin = np.percentile(doubleMutantSignals, 1)
    vmax = np.percentile(doubleMutantSignals, 99)
    vlim = max(abs(vmin), abs(vmax))
    vmin, vmax = -vlim, vlim
    
    # Make heatmap
    plt.figure(figsize=(12,10))
    ax, cax = plotlib.doubleMutant(data, mutantRefVariant, libSeq, 
                                   startPos=startPos, refSignal=targetRefSignal, 
                                   normToRefSignal=norm, coop=coop,
                                   vmin=vmin, vmax=vmax, cbarLabel=cbarLabel,
                                   triangle='lower', invertY=False, linewidth=3)
    cax.tick_params(labelsize=20)
    cax.yaxis.label.set_fontsize(20)
    plotlib.setproperties(title=title, labelfontsize=20, titlefontsize=20)
    
    # Save to file
    if coop:
        figureOutput = 'coop'
    else:
        figureOutput = 'doubleMutant'
    if norm:
        figureNorm = 'norm_'
    else:
        figureNorm = ''
    plt.savefig(figureDir+'/'+names[targetSM]+'_'+names[mutantSM].lower()+'0_'+names[normSM].lower()+'0_'+figureNorm+figureOutput+'.png')
    
    return ax, cax
