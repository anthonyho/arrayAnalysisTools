# Anthony Ho, ahho@stanford.edu, 2/15/2015
# Last update 2/15/2016
"""Python module containing analysis functions for chemical nose project"""


import numpy as np
import pandas as pd
import plotlib


dataDir = 'processedData/'

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

