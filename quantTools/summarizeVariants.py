#!/usr/bin/env python
#
# Summarize a CPseries file by variants given a CPannot file
#
# Anthony Ho, ahho@stanford.edu, 11/29/2016
# Last update 11/29/2016
#
# Work in progress:
# - add different summary metrics


# Import libraries
import os, sys
import argparse
import pandas as pd
from fittinglibs import fileio


def main():

    # Get options and arguments from command line
    parser = argparse.ArgumentParser(description="summarize a CPseries file by variants given a CPannot file")
    parser.add_argument('--outSuffix', help="path to the output file; default is _processed.CPvariant")
    parser.add_argument('annotFilePath', help="path to the CPannot.pkl file")
    parser.add_argument('signalsFilePaths', nargs='*', help="paths to the CPseries to be combined")
    args = parser.parse_args()
    
    # Read CPannot files
    annot = pd.read_pickle(args.annotFilePath)

    # Read signals from CPseries and combine tiles
    allFiles = []
    for i, filePath in enumerate(args.signalsFilePaths):
        allFiles.append(fileio.loadFile(filePath))
    signals = pd.concat(allFiles, axis=1, join='outer')
    signals.columns = range(0, len(allFiles))

    # Define output file path
    if isinstance(args.signalsFilePaths, list):
        (signalDir, signalFilename) = os.path.split(args.signalsFilePaths[0])
    else:
        (signalDir, signalFilename) = os.path.split(args.signalsFilePaths)
    if args.outSuffix is None:
        outputFilePath = os.path.join(signalDir, signalFilename[0:9]+'_processed.CPvariant')
    else:
        outputFilePath = os.path.join(signalDir, signalFilename[0:9]+args.outSuffix)

    # Join signals with CPannot
    signals_with_variants = signals.join(annot, how='inner')

    # Group by and summarize
    signals_by_variants = signals_with_variants.groupby('variant_number').median()
    
    # Write to output file
    signals_by_variants.to_csv(outputFilePath, sep='\t')

    return 1

if __name__ == "__main__":
    main()
