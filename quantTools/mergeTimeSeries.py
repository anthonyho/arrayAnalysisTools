#!/usr/bin/env python
#
# Merge time series data of a tile into a single file
#
# Anthony Ho, ahho@stanford.edu, 1/22/2015
# Last update 1/26/2015


# Import libraries
import os, sys
import argparse
import glob
import pandas as pd
from numpy import pi
import numpy as np
import parselib


def main():

    # Get options and arguments from command line
    parser = argparse.ArgumentParser(description="merge time series data of a tile into a single file")
    parser.add_argument('-r', '--refTime', help="start time of the experiment. If not provided, the start time will be set to the first time point")
    parser.add_argument('numTimepoints', help="number of timepoints")
    parser.add_argument('CPfluorFilePattern', help="pattern of the path to the input CPfluor files. Replace timepoint number with {} and timestamp with *")
    parser.add_argument('outputFilePathNoExt', help="basename of the output files with no extensions")
    args = parser.parse_args()

    # Initialization 

    twoPi = 2 * pi
    
    numTimepoints = int(args.numTimepoints)
    replaceToTimepoint = '##'

    clusterIDColLabels = ['instrumentID', 'runID', 'flowcellID', 'flowcellSide', 'tileID', 'x-pos', 'y-pos']
    fitResultsColLabels = ['fitted', 'amp', 'sigma', 'fittedX', 'fittedY']
    CPfluorColLabels = clusterIDColLabels + fitResultsColLabels

    CPsignals = pd.DataFrame()
    CPsigmas = pd.DataFrame()
    CPtimes = pd.Series()
    
    CPsignalTime = pd.DataFrame(columns=['clusterID','signals','times'])
    CPsigmaTime = pd.DataFrame(columns=['clusterID','sigmas','times'])

    refTime = pd.to_datetime(args.refTime, format='%Y.%m.%d-%H.%M.%S.%f')

    # Output files
    CPsignalTimeFilePath = args.outputFilePathNoExt+".CPsignalTime"
    CPsigmaTimeFilePath = args.outputFilePathNoExt+".CPsigmaTime"


    # Go through all time points
    for timepoint in range(1,numTimepoints+1):
        
        # Get the path to the current CPfluor file
        currTimepoint = str(timepoint)
        currCPfluorFilePattern = args.CPfluorFilePattern.replace(replaceToTimepoint, currTimepoint)
        listOfCurrCPfluorMatchingPattern = glob.glob(currCPfluorFilePattern)

        if len(listOfCurrCPfluorMatchingPattern) > 1:
            # Check if more than one CPfluor file matching pattern exists
            # If so, quit
            print "More than one CPfluor file found in current timepoint!"
            print currCPfluorFilePattern
            print "Quitting now..."
            return 0
        elif len(listOfCurrCPfluorMatchingPattern) == 1:
            # If a unique CPfluor file exists for the current timepoint, proceed
            
            # Get the path to the current CPfluor file
            currCPfluorFile = listOfCurrCPfluorMatchingPattern[0]
            # Load CPfluor file
            currCPfluorData = pd.read_csv(currCPfluorFile, sep=':', names=CPfluorColLabels)
            # Add data from current timepoint to the dataframes
            # If cluster is not fitted, replace value with NaN
            CPsignals[currTimepoint] = ((twoPi * currCPfluorData.amp * currCPfluorData.sigma**2) 
                                        / currCPfluorData.fitted).replace([np.inf, -np.inf], np.nan)
            CPsigmas[currTimepoint] = (currCPfluorData.sigma / currCPfluorData.fitted).replace([np.inf, -np.inf], np.nan)
            # Parse timestamp from CPfluor filename
            CPtimes[currTimepoint] = parselib.parseTimeFromFilename(currCPfluorFile)


    # Compute timepoints relative to the reference timepoint and convert to float64
    if args.refTime == None:
        CPtimes = (CPtimes - CPtimes.iat[0]).apply(lambda x: x / np.timedelta64(1, 's'))
    else:
        CPtimes = (CPtimes - refTime).apply(lambda x: x / np.timedelta64(1, 's'))

    # Get list of columns/indices labels
    listCPsignalsColLabels = list(CPsignals.columns.values)
    listCPsigmasColLabels = list(CPsigmas.columns.values)
    listCPtimesIndLabels = list(CPtimes.index.values)

    # Make CPsignalTime dataframe
    CPsignalTime['clusterID'] = parselib.concatDFColumnsIntoSeries(currCPfluorData, clusterIDColLabels, ':')
    CPsignalTime['signals'] = parselib.concatDFColumnsIntoSeries(CPsignals, listCPsignalsColLabels, ':')
    CPtimesInStr = parselib.concatSeriesIntoString(CPtimes, listCPtimesIndLabels, ':')
    CPsignalTime['times'] = [ CPtimesInStr ]*len(CPsignalTime.index)

    # Make CPsigmaTime dataframe
    CPsigmaTime['clusterID'] = CPsignalTime['clusterID']
    CPsigmaTime['sigmas'] = parselib.concatDFColumnsIntoSeries(CPsigmas, listCPsigmasColLabels, ':')
    CPsigmaTime['times'] = CPsignalTime['times']

    # Write dataframes to files
    CPsignalTime.to_csv(CPsignalTimeFilePath, sep='\t', index=False)
    CPsigmaTime.to_csv(CPsigmaTimeFilePath, sep='\t', index=False)

    return 1 

if __name__ == "__main__":
    main()
