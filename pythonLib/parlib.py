# Anthony Ho, ahho@stanford.edu, 8/4/2016
# Last update 8/9/2016
"""Library functions for embarrassingly parallelization of large lists, Pandas series and dataframes"""


import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import multiprocessing



#
def _splitIntoChunks(listToSplit, numCores):
    lenList = len(listToSplit)
    chunkSize = max(1, lenList // numCores)
    return [listToSplit[i:i + chunkSize] for i in range(0, lenList, chunkSize)]


#
def _funcChunkSeries(dataChunk, func):
    return pd.Series([func(row) for row in dataChunk], index=dataChunk.index)
    

# 
def parallelApply(data, func, numCores, verbose=0, *args, **kwargs):
    '''Parallelizes apply to a large list, Pandas series, or dataframe'''
    if isinstance(data, list):
        results = None
    elif isinstance(data, pd.Series):
        listIndicesChunks = _splitIntoChunks(data.index, numCores)
        listDataChunks = [data.loc[indicesChunk] for indicesChunk in listIndicesChunks]
        listResultsChunks = Parallel(n_jobs=numCores, verbose=verbose)(delayed(_funcChunkSeries)(dataChunk, func) for dataChunk in listDataChunks)
        results = pd.concat(listResultsChunks)
    elif isinstance(data, pd.DataFrame):
        results = None
    
    return results
