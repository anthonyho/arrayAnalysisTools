# Anthony Ho, ahho@stanford.edu, 2/15/2015
# Last update 2/15/2016
"""Python module containing analysis functions for chemical nose project"""


import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import multiprocessing


# 
dataDir = 'processedData/'


# Bootstrap error 
def bootstrap(data, numSamples, statistic, alpha):
    '''Returns bootstrap estimate of 100.0 * (1 - alpha) CI for statistic'''
    n = len(data)
    idx = np.random.randint(0, n, (numSamples, n))
    data_npArray = np.array(data)
    samples = data_npArray[idx]
    stat = np.sort(statistic(samples, 1))
    return (stat[int((alpha / 2.0) * numSamples)],
            stat[int((1 - alpha / 2.0) * numSamples)])


# 
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


#
def normalizeGroupByBootstrap(data_r7, data_QO, data_sm, numClusterCutoff):
    data = pd.DataFrame({'clusterID': data_r7['clusterID'], 
                         'barcode': data_r7['barcode'], 
                         'r7_signals': data_r7['signals'], 
                         'QO_signals': data_QO['signals'], 
                         'sm_signals': data_sm['signals'], 
                         'QO_norm_signals': data_QO['signals'] / data_r7['signals'], 
                         'sm_norm_signals': data_sm['signals'] / data_r7['signals'], 
                         'sm_norm_diff_signals_1': (data_sm['signals'] - data_QO['signals']) / data_r7['signals'], 
                         'sm_norm_diff_signals_2': (data_sm['signals'] - data_QO['signals']) / (data_r7['signals'] - data_QO['signals'])}) 
    data_grouped = data.groupby('barcode', sort=True)
    data_agg = data_grouped['r7_signals', 'QO_signals', 'sm_signals', 'QO_norm_signals', 'sm_norm_signals', 'sm_norm_diff_signals_1'].agg(['size', 'median', 'sem'])
    


    data_agg = data.groupby('barcode', sort=True).agg(['size', 'count', 'median', 'mean', 'sem', 'min', 'max'])
    data_aggPF = data_agg[data_agg['sm_norm_diff_signals_1']['size'] > numClusterCutoff]
    return data, data_agg, data_aggPF

