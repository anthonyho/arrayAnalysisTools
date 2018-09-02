#!/usr/bin/env python
#
# Merge a series of CPfluor files into a single CPseries file
#
# Anthony Ho, ahho@stanford.edu, 7/8/2016
# Last update 9/2/2018


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

    listTimePoint = [str(timepoint) for timepoint in range(1, numTimepoints+1)]

    CPsignals = pd.DataFrame(columns=['clusterID'] + listTimePoint)

    # Output files
    CPsignalsFilePath = args.outputFilePathNoExt+".CPseries"

    # Go through all time points
    for i, currTimepoint in enumerate(listTimePoint):
        
        # Get the path to the current CPfluor file
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

            # Assign clusterID
            if i == 0:
                CPsignals['clusterID'] = parselib.concatDFColumnsIntoSeries(currCPfluorData, clusterIDColLabels, ':')

            # Add data from current timepoint to the dataframes
            # If cluster is not fitted, replace value with NaN
            CPsignals[currTimepoint] = ((twoPi * currCPfluorData.amp * currCPfluorData.sigma**2)
                                        / currCPfluorData.fitted).replace([np.inf, -np.inf], np.nan)
        else:
            CPsignals[currTimepoint] = np.nan
            
        del currCPfluorData

    # Write dataframes to files
    CPsignals.to_csv(CPsignalsFilePath, sep='\t', index=False)

    return 1

if __name__ == "__main__":
    main()
