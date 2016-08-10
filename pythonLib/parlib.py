# Anthony Ho, ahho@stanford.edu, 8/4/2016
# Last update 8/9/2016
"""Library functions for embarrassingly parallelization of large lists, Pandas series and dataframes (by row)"""


import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import multiprocessing



#
def _splitIntoChunks(listToSplit, numChunks):
    lenList = len(listToSplit)
    chunkSize = max(1, lenList // numChunks)
    return [listToSplit[i:i + chunkSize] for i in range(0, lenList, chunkSize)]


#
def _funcChunkSeries(dataChunk, func, *args, **kwargs):
    return dataChunk.apply(func)


#
def _funcChunkDf(dataChunk, func, *args, **kwargs):
    return dataChunk.apply(func, axis=1)


# 
def parallelApply(data, func, numCores, verbose=0, *args, **kwargs):
    '''Parallelizes apply to a large list, Pandas series, or dataframe (row)'''
    if isinstance(data, list):
        results = None
    elif isinstance(data, pd.Series):
        listIndicesChunks = _splitIntoChunks(data.index, numCores)
        listDataChunks = [data.loc[indicesChunk] for indicesChunk in listIndicesChunks]
        listResultsChunks = Parallel(n_jobs=numCores, verbose=verbose)(delayed(_funcChunkSeries)(dataChunk, func, *args, **kwargs) for dataChunk in listDataChunks)
        results = pd.concat(listResultsChunks)
    elif isinstance(data, pd.DataFrame):
        listIndicesChunks = _splitIntoChunks(data.index, numCores)
        listDataChunks = [data.loc[indicesChunk] for indicesChunk in listIndicesChunks]
        listResultsChunks = Parallel(n_jobs=numCores, verbose=verbose)(delayed(_funcChunkDf)(dataChunk, func, *args, **kwargs) for dataChunk in listDataChunks)
        results = pd.concat(listResultsChunks)

    return results
