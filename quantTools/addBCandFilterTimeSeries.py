#!/usr/bin/env python
#
# Add barcodes to time series files, and filter the 
# time series based on filter tag and number of fitted timepoints
#
# Anthony Ho, ahho@stanford.edu, 1/27/2015
# Last update 1/27/2015


# Import libraries
import os, sys
import argparse
import pandas as pd
import numpy as np


# Print analysis summaries
def printSummary(text, number, totalNumber):
    print text, str(number), "("+str(round(float(number)/totalNumber*100, 2))+"%)"
    return


def main():

    # Get options and arguments from command line
    parser = argparse.ArgumentParser(description="add barcodes to time series files, and filter the time series based on filter tag and number of fitted timepoints")
    readGroup = parser.add_mutually_exclusive_group(required=True)
    readGroup.add_argument('-r1', action='store_true', help="read barcode from read 1 if flagged. One of -r1 or -r2 must be present")
    readGroup.add_argument('-r2', action='store_true', help="read barcode from read 2 if flagged. One of -r1 or -r2 must be present")
    parser.add_argument('BCpos', help="position of barcode along read. Use Python notation (e.g. 0:16)")
    parser.add_argument('filterTags', help="filter tags of the clusters to be retained as tagged in filtered tile file. Separate tags with ':'")
    parser.add_argument('maxUnfittedTimepoints', help="maximun number of not fitted timepoints required to retain cluster")
    parser.add_argument('filteredTilePath', help="path to the filtered tile file")
    parser.add_argument('CPtimeSeriesPath', help="path to the CPtimeSeries file")
    parser.add_argument('outputFilePath', help="path to the output file")
    args = parser.parse_args()


    # Define which column to read the barcode from depending on read1 or read2
    if args.r1:
        filteredTileColLabels = ['clusterID', 'filterTags', 'r1', 'q1']
        BCcolLabel = 'r1'
    else:
        filteredTileColLabels = ['clusterID', 'filterTags', 'r1', 'q1', 'r2', 'q2']
        BCcolLabel = 'r2'

    # Load the filtered tile file and the CPtimeSeries file
    filteredTile = pd.read_csv(args.filteredTilePath, sep='\t', usecols=np.arange(len(filteredTileColLabels)), names=filteredTileColLabels)
    CPtimeSeries = pd.read_csv(args.CPtimeSeriesPath, sep='\t')
    signalColIndex = CPtimeSeries.columns[1]

    # Sanity check to make sure filtered tile file and CPtimeSeries file are referring to the same thing
    if len(filteredTile.index) != len(CPtimeSeries.index):
        print "Filtered tile file and CPtimeSeries files have different lengths!"
        print "Quitting..."
        return 0

    # Add barcode to the CPtimeSeries
    startPos = int(args.BCpos.split(':')[0])
    endPos = int(args.BCpos.split(':')[1])
    CPtimeSeries.insert(1, 'barcode', filteredTile[BCcolLabel].str[startPos:endPos])

    
    # Create the regular expression string for searching for filter tags
    listOfFilterTags = args.filterTags.split(':')
    listOfRegexFilterTags = ['\A'+t+'\Z'+'|'+'\A'+t+':'+'|'+':'+t+':'+'|'+':'+t+'\Z' for t in listOfFilterTags]
    regexFilterTags = '|'.join(listOfRegexFilterTags)

    # Create the boolean lists indicating if the cluster passes filtering
    # If the cluster's filter tags contain the right filter tag...
    boolistWithFilterTag = filteredTile['filterTags'].str.contains(regexFilterTags)
    # If the cluster has enough timepoints fitted...
    boolistEnoughFitted = CPtimeSeries[signalColIndex].str.count('nan') <= int(args.maxUnfittedTimepoints)
    # Combine the boolean lists...
    boolistPassingFilter = boolistWithFilterTag & boolistEnoughFitted 


    # Print summaries
    print "Analyzing the following files:"
    print "  "+args.filteredTilePath
    print "  "+args.CPtimeSeriesPath
    printSummary("  Total number of clusters analyzed:", len(CPtimeSeries.index), len(CPtimeSeries.index))
    printSummary("  Number of clusters with the correct filter tags:", np.sum(boolistWithFilterTag), len(CPtimeSeries.index))
    printSummary("  Number of clusters with enough timepoints fitted:", np.sum(boolistEnoughFitted), len(CPtimeSeries.index))
    printSummary("  Number of clusters with the correct filter tags and enough timepoints fitted:", 
                 np.sum(boolistPassingFilter), len(CPtimeSeries.index))

    # Write the doubly filtered time series with barcodes to output file
    CPtimeSeries[boolistPassingFilter].to_csv(args.outputFilePath, sep='\t', index=False)
    
    return 1 

if __name__ == "__main__":
    main()
