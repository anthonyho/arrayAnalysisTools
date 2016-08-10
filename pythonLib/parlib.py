# Anthony Ho, ahho@stanford.edu, 8/4/2016
# Last update 8/9/2016
"""Library functions for embarrassingly parallelization applying functions to large lists, Pandas series, df (by row), and grouped df"""


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
def _funcChunkList(dataChunk, func, *args, **kwargs):
    return [func(item, *args, **kwargs) for item in dataChunk]


#
def _funcChunkSeries(dataChunk, func, *args, **kwargs):
    return dataChunk.apply(func, args=args, **kwargs)


#
def _funcChunkDf(dataChunk, func, *args, **kwargs):
    return dataChunk.apply(func, axis=1, args=args, **kwargs)


#
def _funcChunkGroupApply(dataChunk, func, indexKeys, *args, **kwargs):

#    listNames = [name for name, group in dataChunk]
    listResults = [func(group, *args, **kwargs) for name, group in dataChunk]
    results = pd.concat(listResults)

#    if isinstance(listNames[0], tuple):
#        results = pd.concat(listResults)#, index=pd.MultiIndex.from_tuples(listNames, names=indexKeys))
#    else:
#        results = pd.concat(listResults)#, index=listNames)
#        results.index.name = indexKeys

    return results


#
def _funcChunkGroupAgg(dataChunk, func, indexKeys, *args, **kwargs):

    listNames = [name for name, group in dataChunk]
    listResults = [func(group, *args, **kwargs) for name, group in dataChunk]

    if isinstance(listNames[0], tuple):
        results = pd.DataFrame(listResults, index=pd.MultiIndex.from_tuples(listNames, names=indexKeys))
    else:
        results = pd.DataFrame(listResults, index=listNames)
        results.index.name = indexKeys

    return results


#
def _funcChunkGroupAggCompact(dataChunk, func, indexKeys, *args, **kwargs):

    listNames = [name for name, group in dataChunk]
    listResults = [func(group, *args, **kwargs) for name, group in dataChunk]

    if isinstance(listNames[0], tuple):
        results = pd.Series(listResults, index=pd.MultiIndex.from_tuples(listNames, names=indexKeys))
    else:
        results = pd.Series(listResults, index=listNames)
        results.index.name = indexKeys

    return results


# 
def parallelApply(data, func, numCores, verbose=0, *args, **kwargs):
    '''Apply a function in parallel to a large list, Pandas series, df (by row), or grouped df'''

    # List
    if isinstance(data, list):
        listDataChunks = _splitIntoChunks(data, numCores)
        listResultsChunks = Parallel(n_jobs=numCores, verbose=verbose)(delayed(_funcChunkList)(dataChunk, func, *args, **kwargs) 
                                                                       for dataChunk in listDataChunks)
        results = [item for sublist in listResultsChunks for item in sublist]

    # Pandas series
    elif isinstance(data, pd.Series):
        listIndicesChunks = _splitIntoChunks(data.index, numCores)
        listDataChunks = [data.loc[indicesChunk] for indicesChunk in listIndicesChunks]
        listResultsChunks = Parallel(n_jobs=numCores, verbose=verbose)(delayed(_funcChunkSeries)(dataChunk, func, *args, **kwargs) 
                                                                       for dataChunk in listDataChunks)
        results = pd.concat(listResultsChunks)

    # Pandas dataframe
    elif isinstance(data, pd.DataFrame):
        listIndicesChunks = _splitIntoChunks(data.index, numCores)
        listDataChunks = [data.loc[indicesChunk] for indicesChunk in listIndicesChunks]
        listResultsChunks = Parallel(n_jobs=numCores, verbose=verbose)(delayed(_funcChunkDf)(dataChunk, func, *args, **kwargs) 
                                                                       for dataChunk in listDataChunks)
        results = pd.concat(listResultsChunks)

    # Pandas grouped dataframe
    # Aggregate only for now, compatible with groupby multiple columns
    elif isinstance(data, pd.core.groupby.DataFrameGroupBy):
        # Get name of the columns being grouped
        if (isinstance(data.keys, list) or isinstance(data.keys, tuple)) and len(data.keys) == 1:
            indexKeys = data.keys[0]
        else:
            indexKeys = data.keys
        # Convert groupby object into list of tuples of names and groups
        listNamesGroups = list(data)
        # Test on the first group to see the type of output from func(grouped)
        returnType = func(listNamesGroups[0][1], *args, **kwargs)
        # Break into chunks, apply func, concat results
        listDataChunks = _splitIntoChunks(listNamesGroups, numCores)
        if isinstance(returnType, pd.DataFrame):
            listResultsChunks = Parallel(n_jobs=numCores, verbose=verbose)(delayed(_funcChunkGroupApply)(dataChunk, func, indexKeys, *args, **kwargs) 
                                                                           for dataChunk in listDataChunks)
        elif isinstance(returnType, pd.Series):
            listResultsChunks = Parallel(n_jobs=numCores, verbose=verbose)(delayed(_funcChunkGroupAgg)(dataChunk, func, indexKeys, *args, **kwargs) 
                                                                           for dataChunk in listDataChunks)
        else:
            listResultsChunks = Parallel(n_jobs=numCores, verbose=verbose)(delayed(_funcChunkGroupAggCompact)(dataChunk, func, indexKeys, *args, **kwargs) 
                                                                           for dataChunk in listDataChunks)
        results = pd.concat(listResultsChunks)

    return results


# ??
# determine series vs df output
