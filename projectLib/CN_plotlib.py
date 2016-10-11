# Anthony Ho, ahho@stanford.edu, 9/1/2015
# Last update 9/1/2016
"""Python module containing plot functions for chemical nose project"""


import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import plotlib
import liblib
from fittinglibs import plotting, processresults
import CN_analysislib


# Plot fitted binding curve of a single variant
# Borrowed from Sarah's pipeline with modifications to show Kd instead of dG, 
# and allow for using current axis
def plotBindingCurve(affinityData, variant, annotate=True, ax=None):
    '''Plot a binding curve of a particular variant'''
    # Load data
    subSeries = affinityData.getVariantBindingSeries(variant)
    if len(subSeries)==0:
        print 'No fluorescence data associated with variant %s'%str(variant)
        return
    concentrations = affinityData.x
    variant_table = affinityData.variant_table
    
    # Plot
    if ax is None:
        fig = plt.figure(figsize=(4, 4))
        ax = fig.add_subplot(111)
    plotting.plotFitCurve(concentrations, subSeries, variant_table.loc[variant], ax=ax)

    plotlib.setproperties(xlabel='Concentration (nM)', ylabel='Fluorescence', 
                          labelfontsize=16, tickfontsize=16)

    # Annonate
    if annotate:
        names = ['dG', 'dG_lb', 'dG_ub', 'fmax', 'fmax_lb', 'fmax_ub', 'numTests', 'pvalue', 'rsq', 'flag']
        vec = pd.Series([processresults.getValueInTable(variant_table.loc[variant], name) for name in names], index=names)
            #vec.fillna('', inplace=True)
        annotationText = ['Kd = {:4.2f} uM ({:4.2f}, {:4.2f})'.format(liblib.dGtoKd(vec.dG, 'uM'),
                                                                   liblib.dGtoKd(vec.dG_lb, 'uM'),
                                                                   liblib.dGtoKd(vec.dG_ub, 'uM')),
                          'fmax = {:4.2f} ({:4.2f}, {:4.2f})'.format(vec.fmax,
                                                                     vec.fmax_lb,
                                                                     vec.fmax_ub),
                          'Nclusters = {:4.0f}'.format(vec.numTests),
                          'pvalue = {:.1e}'.format(vec.pvalue),
                          'average Rsq = {:4.2f}'.format(vec.rsq),
                          'fmax enforced = {}'.format(vec.flag)
                          ]
        ax.annotate('\n'.join(annotationText), xy=(0.05, .95), xycoords='axes fraction',
                    horizontalalignment='left', verticalalignment='top', fontsize=11)

    return ax


# Plot fitted binding curves of a single variant across different targets
# Borrowed from Sarah's pipeline with modifications to show Kd instead of dG, 
# and allow for using current axis
def plotBindingCurvesAcrossTargets(affinityDataList, variant, 
                                  colors=None, names=None, norm=False, ax=None):
    '''Plot a binding curve of a particular variant'''

    if colors is None:
        colors = sns.color_palette("Paired", 12)

    if ax is None:
        fig = plt.figure(figsize=(4, 4))
        ax = fig.add_subplot(111)

    for i, affinityData in enumerate(affinityDataList):
        
        # Load data
        subSeries = affinityData.getVariantBindingSeries(variant)
        if len(subSeries)==0:
            print 'No fluorescence data associated with variant %s'%str(variant)
            return
        concentrations = affinityData.x
        variant_table = affinityData.variant_table
        variant_data = variant_table.loc[variant]

        if norm:
            fmax = variant_data['fmax']
            subSeries = subSeries.copy() / fmax
            variant_data = variant_data.copy()
            for field in ['fmax_init', 'fmin_init', 
                          'fmax_lb', 'fmax', 'fmax_ub', 
                          'fmin_lb', 'fmin', 'fmin_ub']:
                variant_data[field] = variant_data[field] / fmax
        
        # Plot
        plotting.plotFitCurve(concentrations, subSeries, variant_data, ax=ax)
        
        ax.lines[i*4+3].set_color(colors[i])
        ax.lines[i*4+3].set_linewidth(2)
        if names is not None:
            ax.lines[i*4+3].set_label(names[i])
    
    if norm:
        ylabel = 'Normalized fluorescence'
    else:
        ylabel = 'Fluorescence'
    plotlib.setproperties(xlabel='Concentration (nM)', ylabel=ylabel, 
                          labelfontsize=16, tickfontsize=16, 
                          legend=names, legendloc=2, legendfontsize=12)

    return ax


# Plot fitted binding curves of different variant on the same target
# Borrowed from Sarah's pipeline with modifications to show Kd instead of dG, 
# and allow for using current axis
def plotBindingCurvesAcrossVariants(affinityData, variantList, 
                                    colors=None, norm=False, ax=None):
    '''Plot a binding curve of a particular variant'''

    if colors is None:
        colors = sns.color_palette("Paired", 12)

    if ax is None:
        fig = plt.figure(figsize=(4, 4))
        ax = fig.add_subplot(111)

    for i, variant in enumerate(variantList):
        
        # Load data
        subSeries = affinityData.getVariantBindingSeries(variant)
        if len(subSeries)==0:
            print 'No fluorescence data associated with variant %s'%str(variant)
            return
        concentrations = affinityData.x
        variant_table = affinityData.variant_table
        variant_data = variant_table.loc[variant]
        
        if norm:
            fmax = variant_data['fmax']
            subSeries = subSeries.copy() / fmax
            variant_data = variant_data.copy()
            for field in ['fmax_init', 'fmin_init', 
                          'fmax_lb', 'fmax', 'fmax_ub', 
                          'fmin_lb', 'fmin', 'fmin_ub']:
                variant_data[field] = variant_data[field] / fmax
        
        # Plot
        plotting.plotFitCurve(concentrations, subSeries, variant_data, ax=ax)
        
        ax.lines[i*4+3].set_color(colors[i])
        ax.lines[i*4+3].set_linewidth(2)

    if norm:
        ylabel = 'Normalized fluorescence'
    else:
        ylabel = 'Fluorescence'    
    plotlib.setproperties(xlabel='Concentration (nM)', ylabel=ylabel, 
                          labelfontsize=16, tickfontsize=16)

    return ax


# Plot double mutant/cooperativity heatmap for one small molecule
def _doubleMutant(variant_table_list, mutantRefVariant, mutantSM, targetSM, normSM, names, norm=False, coop=False):
    # Define variables
    data = variant_table_list[targetSM]['dG']
    
    #####
    targetRefVariant = 'gtGGATTTTCCgcatACGAAGTtgTCCCGA'.upper()
    normRefVariant = 'gtGGATTTTCCgcatACGAAGTtgTCCCGA'.upper()
    #####

#    mutantRefVariant = variant_table_list[mutantSM]['dG'].idxmin()
#    targetRefVariant = variant_table_list[targetSM]['dG'].idxmin()
#    normRefVariant = variant_table_list[normSM]['dG'].idxmin()
    if norm:
        targetRefSignal = variant_table_list[targetSM]['dG'][targetRefVariant]
        normRefSignal = variant_table_list[normSM]['dG'][normRefVariant]
    else:
        targetRefSignal = 0 ####
        normRefSignal = 0   ####
    refLibSeq = CN_analysislib.returnRefSeq(mutantRefVariant)
    startPos = 14 #####
    title = 'Mutants of '+names[mutantSM].lower()+'0, in the presence of '+names[targetSM]
    if norm and coop:
        cbarLabel = 'Cooperativity normalized to '+names[targetSM].lower()+'0'
    elif norm and not coop:
        cbarLabel = 'dG normalized to '+names[targetSM].lower()+'0'
    elif not norm and coop:
        cbarLabel = 'Cooperativity'
    else:
        cbarLabel = 'dG'
    
    # Compute vmin and vmax manually (assuming robust=True) so that the color bar
    # scale is the same across the same dataRefVariant
    #normRefData = variant_table_list[normSM]['dG'] / normRefSignal
    normRefData = variant_table_list[normSM]['dG'] - normRefSignal
    doubleMutantSignals, mutantLabels = plotlib.doubleMutantMatrix(normRefData, normRefVariant, refLibSeq, startPos, coop)
    doubleMutantSignals = doubleMutantSignals[~np.isnan(doubleMutantSignals)]
    vmin = np.percentile(doubleMutantSignals, 1) ####
    vmax = np.percentile(doubleMutantSignals, 99) ####
    vlim = max(abs(vmin), abs(vmax))
    vmin, vmax = -vlim, vlim
    
    # Make heatmap
    plt.figure(figsize=(12,10))
    ax, cax = plotlib.doubleMutant(data, mutantRefVariant, refLibSeq, 
                                   startPos=startPos, refSignal=targetRefSignal, 
                                   normToRefSignal=norm, coop=coop,
###                                   vmin=vmin, vmax=vmax, cbarLabel=cbarLabel,
                                   cbarLabel=cbarLabel, center=1,
                                   triangle='lower', invertY=False, linewidth=3, cmap="RdYlBu") ###
    cax.tick_params(labelsize=20)
    cax.yaxis.label.set_fontsize(20)
    plotlib.setproperties(title=title, labelfontsize=20, titlefontsize=20, 
                          xticklabelrot=90, yticklabelrot=0)
    
    # Save to file
#    if coop:
#        figureOutput = 'coop'
#    else:
#        figureOutput = 'doubleMutant'
#    if norm:
#        figureNorm = 'norm_'
#    else:
#        figureNorm = ''
#    plt.savefig(figureDir+'/'+names[targetSM]+'_'+names[mutantSM].lower()+'0_'+names[normSM].lower()+'0_'+figureNorm+figureOutput+'.png')
    
    return ax, cax
