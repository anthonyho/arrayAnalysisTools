# Anthony Ho, ahho@stanford.edu, 2/3/2015
# Last update 6/12/2015
"""Library of parsing tools specific to the array analysis pipeline"""


import os
import numpy as np
import pandas as pd


# Split a column containing separator-separated values into individual
# columns and assign to new numpy array
def splitConcatedDFColumnIntoNDarray(column, separator=':'):
    return np.array(column.str.split(separator).tolist()).astype(float)


# Split a column containing separator-separated values into individual
# columns and assign to new dataframe
def splitConcatedDFColumnIntoDF(column, separator=':', string=False):
    if string:
        return pd.DataFrame(column.str.split(separator).tolist())
    else:
        return pd.DataFrame(column.str.split(separator).tolist()).convert_objects(convert_numeric=True)


# Concatenate designated columns of a datadrame into a series of
# strings, separated by separator
def concatDFColumnsIntoSeries(df, columnsLabels, separator=':'):
    series = pd.Series(df[columnsLabels[0]].map(str))
    for currCol in columnsLabels[1:]:
        series = series+separator+df[currCol].map(str)
    return series


# Concatenate designated elements of a series into a string
# separated by separator
def concatSeriesIntoString(series, indexLabels, separator=':'):
    string = str(series[indexLabels[0]])
    for currElement in indexLabels[1:]:
        string = string+separator+str(series[currElement])
    return string


# Parse timestamp from the last part of the filename as datetime64[ns] objects
def parseTimeFromFilename(fileFullPath):
    (_, filename) = os.path.split(fileFullPath)
    (fileBasename, _) = os.path.splitext(filename)
    timestampStr = fileBasename.split('_')[-1]
    timestamp = pd.to_datetime(timestampStr, format='%Y.%m.%d-%H.%M.%S.%f')
    return timestamp


# Get a specific column from each of the df in a list 
# and concat them into a new df
def extractConcatCol(listDf, col, join='outer', names=None):
    listSeries = [df[col] for i, df in enumerate(listDf)]
    concatedCols = pd.concat(listSeries, axis=1, join=join)
    if names is not None:
        concatedCols.columns = names
    return concatedCols
