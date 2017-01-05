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


# Convert delta G in kcal/mol to Kd
# default Kd in uM
# default temperature at 20C
def dGtoKd(data, unit='uM', T=20):
    '''Computes delta G in kcal/mol to Kd '''
    # Compute RT
    RT = 0.0019872036 * (T + 273.15)

    # Set unit multiplier
    if unit == 'pM':
        multiplier = 1e12
    elif unit == 'nM':
        multiplier = 1e9
    elif unit == 'uM':
        multiplier = 1e6
    elif unit == 'mM':
        multiplier = 1e3
    elif unit == 'M':
        multiplier = 1
    else:
        raise ValueError('Unit \"'+unit+'\" not supported!')
    
    try: 
        return np.exp(data.astype(float) / RT) * multiplier
    except AttributeError:
        return np.exp(data / RT) * multiplier

