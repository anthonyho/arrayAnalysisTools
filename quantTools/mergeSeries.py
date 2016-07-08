#!/usr/bin/env python
#
# Merge a series of CPfluor files into a single CPseries file
#
# Anthony Ho, ahho@stanford.edu, 7/8/2016
# Last update 7/8/2016


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
    parser = argparse.ArgumentParser(description="merge a series of CPfluor files into a single CPseries")
    parser.add_argument('-t', '--time', help="if enabled, add an additional column with comma separated time for each CPfluor (default: false)", action='store_true', default=False)
    parser.add_argument('-r', '--refTime', help="start time of the experiment. If not provided, the start time will be set to the first time point")
    parser.add_argument('numTimepoints', help="number of timepoints")
    parser.add_argument('CPfluorFilePattern', help="pattern of the path to the input CPfluor files. Replace timepoint number with ## and timestamp with *")
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
    if args.time:
        CPtimes = pd.Series()
        CPseries = pd.DataFrame(columns=['clusterID', 'signals', 'times'])
        refTime = pd.to_datetime(args.refTime, format='%Y.%m.%d-%H.%M.%S.%f')
    else:
        CPseries = pd.DataFrame(columns=['clusterID', 'signals'])

    # Output files
    CPseriesFilePath = args.outputFilePathNoExt+".CPseries"

    # Go through all time points
    for timepoint in range(1, numTimepoints+1):

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
            # Parse timestamp from CPfluor filename
            if args.time:
                CPtimes[currTimepoint] = parselib.parseTimeFromFilename(currCPfluorFile)

    # Compute timepoints relative to the reference timepoint and convert to float64
    if args.time:
        if args.refTime is None:
            CPtimes = (CPtimes - CPtimes.iat[0]).apply(lambda x: x / np.timedelta64(1, 's'))
        else:
            CPtimes = (CPtimes - refTime).apply(lambda x: x / np.timedelta64(1, 's'))

    # Get list of columns/indices labels
    listCPsignalsColLabels = list(CPsignals.columns.values)
    if args.time:
        listCPtimesIndLabels = list(CPtimes.index.values)

    # Make CPsignalTime dataframe
    CPseries['clusterID'] = parselib.concatDFColumnsIntoSeries(currCPfluorData, clusterIDColLabels, ':')
    CPseries['signals'] = parselib.concatDFColumnsIntoSeries(CPsignals, listCPsignalsColLabels, ':')
    if args.time:
        CPtimesInStr = parselib.concatSeriesIntoString(CPtimes, listCPtimesIndLabels, ':')
        CPseries['times'] = [CPtimesInStr]*len(CPsignalTime.index)

    # Write dataframes to files
    CPseries.to_csv(CPseriesFilePath, sep='\t', index=False)

    return 1

if __name__ == "__main__":
    main()
