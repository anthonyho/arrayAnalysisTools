#!/usr/bin/env python
#
# Compress duplicated sequence variants into unique sequence variants
# count number of duplicates, and compute the corresponding median/mean
# of the various fields
#
# Anthony Ho, ahho@stanford.edu, 2/24/2015
# Last update 3/3/2015


# Import libraries
import os, sys
import argparse
import pandas as pd
import numpy as np


def main():

    # Get options and arguments from command line
    parser = argparse.ArgumentParser(description="compress variants")
    parser.add_argument('-c', '--count', action='store_true', help="set to true if the count of sequences within a barcode block is at the end of annotation (default=false)")
    parser.add_argument('inputFilePath', help="path to the input file")
    parser.add_argument('outputFilePrefix', help="path to the output file without extension")
    args = parser.parse_args()


    # --- Initialize --- #

    # Define a list of prefixes of the column name to compute mean, median and SD
    listColPrefix = ['params', 'paramSEs', 'RSS', 'reChi2', 'SER']

    # Define output file extension
    groupedClustersExt = '.CPfitGrouped'
    unqClustersExt = '.CPfitUnq'


    # --- Read input file and groupby annotations --- #

    # Read clusters from inputFile
    allClusters = pd.read_csv(args.inputFilePath, sep='\t')

    # Drop the last colon-separated field (barcode block count) from annotations if count mode is on
    if args.count:
        allClusters['annotation'] = allClusters['annotation'].replace(":[0-9]+$", ":", regex=True)  # need to check regex

    # Group by annotation
    grouped = allClusters.groupby('annotation')

    # Compute number of count of each unique sequence variant
    groupedCounts = grouped.size()
    groupedCounts.name = 'count'


    # --- Write all clusters to file grouped by annotations and sorted by count --- #

    allClustersIndAnnt = allClusters.set_index('annotation')
    allClustersIndAnnt.insert(1, 'count', groupedCounts)
    allClustersGroupedSorted = allClustersIndAnnt.reset_index().sort(['count', 'annotation'], ascending=[False, True])

    # Write to output file
    allClustersGroupedSorted.to_csv(args.outputFilePrefix+groupedClustersExt, sep='\t', index=False)


    # --- Compute unique sequence statistics and write to file --- #

    # Compute the list of columns to compute mean, median and SD
    listCol = []
    for prefix in listColPrefix:
        listCol.extend([col for col in allClusters.columns.values if col.startswith(prefix)])

    # Compute medians, means, and standard deviation
    groupedMedians = grouped[listCol].agg(np.median).add_suffix('.median')
    groupedMeans = grouped[listCol].agg(np.mean).add_suffix('.mean')
    groupedSDs = grouped[listCol].agg(np.std).add_suffix('.sd')

    # Concat the group counts, medians, means, and std into the same dataframe
    uniqueClusters = pd.concat([groupedCounts, groupedMedians, groupedMeans, groupedSDs], axis=1).reset_index()
    # Sort in descending order by the group count
    uniqueClusters = uniqueClusters.sort('count', ascending=False)

    # Write to output file
    uniqueClusters.to_csv(args.outputFilePrefix+unqClustersExt, sep='\t', index=False)

    return 1

if __name__ == "__main__":
    main()
