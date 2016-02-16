# Anthony Ho, ahho@stanford.edu, 2/15/2015
# Last update 2/15/2016
"""Python module containing analysis functions for chemical nose project"""


import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import multiprocessing
import liblib


# Combine raw data and computed data into a single dataframe
def normalize(data_r7, data_QO, data_sm):
    data = pd.DataFrame({'clusterID': data_r7['clusterID'], 
                         'barcode': data_r7['barcode'], 
                         'r7_signals': data_r7['signals'], 
                         'QO_signals': data_QO['signals'], 
                         'sm_signals': data_sm['signals'], 
                         'QO_norm_signals': data_QO['signals'] / data_r7['signals'], 
                         'sm_norm_signals': data_sm['signals'] / data_r7['signals'], 
                         'sm_norm_diff_signals_1': (data_sm['signals'] - data_QO['signals']) / data_r7['signals'], 
                         'sm_norm_diff_signals_2': (data_sm['signals'] - data_QO['signals']) / (data_r7['signals'] - data_QO['signals'])}) 
    return data


# Group by and summarize data
def summarize(data, numCores):
    # Group by barcodes
    data_grouped = data.groupby('barcode', sort=True)
    # Compute summary statistics for the variables not being bootstrapped
    data_agg_notbs = data_grouped['r7_signals', 'QO_signals', 'sm_signals', 'QO_norm_signals', 'sm_norm_signals', 'sm_norm_diff_signals_1'].agg(['size', 'median', 'sem'])
    # Compute summary statistics for the variables being bootstrapped
    data_agg_bs_other = data_grouped['sm_norm_diff_signals_2'].agg(['size', 'median'])
    data_agg_bs_err = liblib.aggParallel(data_grouped['sm_norm_diff_signals_2'], liblib.bootstrap, 'bserr', numCores, 5, np.median, 0.05)
    data_agg_bs = pd.concat({'sm_norm_diff_signals_2': pd.concat([data_agg_bs_other, data_agg_bs_err], axis=1)}, axis=1)
    # Combine summary statistics 
    data_agg = pd.concat([data_agg_notbs, data_agg_bs], axis=1)

    return data_agg


# Filter data
def filter(data_agg, numClusterCutoff): 
    data_aggPF = data_agg[data_agg['sm_norm_diff_signals_2']['size'] > numClusterCutoff]
    return data_aggPF


# Deprecated - use normalize(), summarize(), and filter() instead
def normalizeGroupBy(data_r7, data_QO, data_sm, numClusterCutoff):
    data = pd.DataFrame({'clusterID': data_r7['clusterID'], 
                         'barcode': data_r7['barcode'], 
                         'r7_signals': data_r7['signals'], 
                         'QO_signals': data_QO['signals'], 
                         'sm_signals': data_sm['signals'], 
                         'QO_norm_signals': data_QO['signals'] / data_r7['signals'], 
                         'sm_norm_signals': data_sm['signals'] / data_r7['signals'], 
                         'sm_norm_diff_signals_1': (data_sm['signals'] - data_QO['signals']) / data_r7['signals'], 
                         'sm_norm_diff_signals_2': (data_sm['signals'] - data_QO['signals']) / (data_r7['signals'] - data_QO['signals'])}) 
    data_agg = data.groupby('barcode', sort=True).agg(['size', 'count', 'median', 'mean', 'sem', 'min', 'max'])
    data_aggPF = data_agg[data_agg['sm_norm_diff_signals_1']['size'] > numClusterCutoff]
    return data, data_agg, data_aggPF


