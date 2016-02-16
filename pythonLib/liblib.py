# Anthony Ho, ahho@stanford.edu, 2/15/2016
# Last update 2/15/2016
"""Miscellaneous library functions"""


import pandas as pd
from joblib import Parallel, delayed
import multiprocessing


# Compute a summary statistic to a grouped Pandas dataframe in parallel
# and return a Pandas series
def aggParallel(dfGrouped, func, name, numCores, *args, **kwargs):
    return pd.Series(Parallel(n_jobs=numCores)(delayed(func)(group, *args, **kwargs) for name, group in dfGrouped), 
                     name=name,
                     index=dfGrouped.indices)
