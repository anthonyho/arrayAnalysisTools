# Anthony Ho, ahho@stanford.edu, 2/3/2015
# Last update 2/9/2015
""" Library of parsing tools specific to the array analysis pipeline """


#import subprocess
import numpy as np
import pandas as pd
#import string


# Concatenate designated columns of a datadrame into a series of 
# strings, separated by separator
def concatDFColumnsIntoSeries(df, columnsLabels, separator):
    s = pd.Series(df[columnsLabels[0]].map(str))
    for currCol in columnsLabels[1:]:
        s = s+separator+df[currCol].map(str)
    return s    

# Concatenate designated elements of a series into a string
# separated by separator
def concatSeriesIntoString(series, indexLabels, separator):
    string = str(series[indexLabels[0]])
    for currElement in indexLabels[1:]:
        string = string+separator+str(series[currElement])
    return string

# Parse timestamp from the last part of the filename as datetime64[ns] objects
def parseTimeFromFilename(fileFullPath):
    (dirPath, filename) = os.path.split(fileFullPath)
    (fileBasename, fileExt) = os.path.splitext(filename)
    timestampStr = fileBasename.split('_')[-1]
    timestamp = pd.to_datetime(timestampStr, format='%Y.%m.%d-%H.%M.%S.%f')
    return timestamp

