# Anthony Ho, ahho@stanford.edu, 2/15/2016
# Last update 8/4/2016
"""Miscellaneous library functions"""


import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import multiprocessing


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


# Compute a summary statistic to a grouped Pandas dataframe in parallel
# and return a Pandas series
# Deprecated - use functions in parlib instead
def aggParallel(dfGrouped, func, name, numCores, *args, **kwargs):
    '''Computes a summary statistics to a grouped Pandas dataframe in parallel'''
    return pd.Series(Parallel(n_jobs=numCores)(delayed(func)(group, *args, **kwargs) for name, group in dfGrouped),
                     name=name,
                     index=sorted(dfGrouped.indices.keys()))
