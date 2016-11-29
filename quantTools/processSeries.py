#!/usr/bin/env python
#
# Combine signals from multiple CPseries into one file, filter by clusters in CPannot, 
# and apply custom function to transform signals
#
# Anthony Ho, ahho@stanford.edu, 8/18/2016
# Last update 8/18/2016


# Import libraries
import os, sys
import argparse
import pandas as pd


def parseFuncFromFile(funcFilePath):
    
    # Import func from funcFilePath
    (funcDir, funcFilename) = os.path.split(funcFilePath)
    (funcFileBasename, _) = os.path.splitext(funcFilename)
    
    sys.path.insert(0, funcDir)
    func_module = __import__(funcFileBasename)

    return vars(func_module)['func']


def main():

    # Get options and arguments from command line
    parser = argparse.ArgumentParser(description="combine, filter, and transform signals from different CPseries into one file")
    parser.add_argument('--outSuffix', help="suffix of the output file; _processed.CPseries.pkl by defualt")
    parser.add_argument('annotFilePath', help="path to the CPannot.pkl file")
    parser.add_argument('funcFilePath', help="path to the file specifying the transform function")
    parser.add_argument('signalsFilePaths', nargs='*', help="paths to the CPseries to be combined")
    args = parser.parse_args()
    
    # Read func and CPannot files
    func = parseFuncFromFile(args.funcFilePath)
    annot = pd.read_pickle(args.annotFilePath)

    # Read signals from CPseries and combine tiles
    allClusters = []
    for i, tilePath in enumerate(args.signalsFilePaths):
        allClusters.append(pd.read_csv(tilePath, sep='\t'))
    signals = pd.concat(allClusters, axis=0, join='outer', ignore_index=True)

    # Define output file path
    if isinstance(args.signalsFilePaths, list):
        (signalDir, signalFilename) = os.path.split(args.signalsFilePaths[0])
    else:
        (signalDir, signalFilename) = os.path.split(args.signalsFilePaths)
    if args.outSuffix is None:
        outputFilePath = os.path.join(signalDir, signalFilename[0:9]+'_processed.CPseries.pkl')
    else:
        outputFilePath = os.path.join(signalDir, signalFilename[0:9]+args.outSuffix)

    # Filter clusters by CPannot
    signals_filtered = signals[signals['clusterID'].isin(annot.index)].set_index('clusterID')
    
    # Apply custom transformation function
    signals_filtered_transformed = func(signals_filtered)
    
    # Write to output file
    signals_filtered_transformed.to_pickle(outputFilePath)

    return 1

if __name__ == "__main__":
    main()
