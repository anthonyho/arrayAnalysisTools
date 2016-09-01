# Anthony Ho, ahho@stanford.edu, 6/4/2016
# Last update 6/4/2016
"""Python module containing analysis functions for chemical nose project"""


import numpy as np
import pandas as pd
import liblib


# Combine raw data and computed data into a single dataframe
def normalize(data):
    data_norm = pd.DataFrame({'clusterID': data['clusterID'], 
                              'tag': data['tag'], 
                              'barcode': data['barcode'], 
                              'seq': data['seq'], 
                              'allCluster': data['allCluster'], 
                              'quenched': data['quenched'],
                              'switch': data['switch'], 
                              'quenched_norm': data['quenched'] / data['allCluster'], 
                              'switch_norm': data['switch'] / data['allCluster'], 
                              'diff': (data['switch'] - data['quenched']) / (data['allCluster'] - data['quenched']), 
                              }) 
    return data_norm


# Group by and summarize data
def summarize(data_norm):
#def summarize(data, numCores, numSamples, alpha):
    # Group by barcodes
    data_norm_grouped = data_norm.groupby('seq', sort=True)
    data_agg = data_norm_grouped['allCluster', 'quenched', 'switch', 
                                 'quenched_norm', 'switch_norm', 'diff'].agg(['size', 'median', 'sem'])

    # Compute summary statistics for the variables not being bootstrapped
#    data_agg_notbs = data_norm_grouped['r7_signals', 'QO_signals', 'sm_signals', 'QO_norm_signals', 'sm_norm_signals', 'sm_norm_diff_signals_1'].agg(['size', 'median', 'sem'])
    # Compute summary statistics for the variables being bootstrapped
#    data_agg_bs_other = data_grouped['sm_norm_diff_signals_2'].agg(['size', 'median'])
#    data_agg_bs_err = liblib.aggParallel(data_grouped['sm_norm_diff_signals_2'], liblib.bootstrap, 'bserr', numCores, numSamples, np.median, alpha)
#    data_agg_bs = pd.concat({'sm_norm_diff_signals_2': pd.concat([data_agg_bs_other, data_agg_bs_err], axis=1)}, axis=1)
    # Combine summary statistics 
#    data_agg = pd.concat([data_agg_notbs, data_agg_bs], axis=1)

    return data_agg


# Filter data
def filter(data_agg, numClusterCutoff): 
    data_aggPF = data_agg[data_agg['switch']['size'] > numClusterCutoff]
    return data_aggPF

