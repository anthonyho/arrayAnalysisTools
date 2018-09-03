#!/usr/bin/env python
#
# Summarize a CPseries file by variants given a CPannot file
#
# Anthony Ho, ahho@stanford.edu, 11/29/2016
# Last update 9/2/2018
#
# Work in progress:
# - add different summary metrics


# Import libraries
import os, sys
import argparse
import pandas as pd


def main():

    # Get options and arguments from command line
    parser = argparse.ArgumentParser(description="load and summarize all CPseries files by variants given a CPannot file")
    parser.add_argument('annotFilePath', help="path to the CPannot.pkl file")
    parser.add_argument('signalsFilePaths', nargs='*', help="paths to the CPseries to be combilned")
    args = parser.parse_args()
    
    # Read CPannot files
    annot = pd.read_pickle(args.annotFilePath)

    # Read signals from CPseries and combine tiles
    allFiles = []
    for i, filePath in enumerate(args.signalsFilePaths):
        allFiles.append(pd.read_csv(filePath, sep='\t'))
    signals = pd.concat(allFiles, axis=0, join='outer').set_index('clusterID')

    # Define output file path
    if isinstance(args.signalsFilePaths, list):
        (signalDir, signalFilename) = os.path.split(args.signalsFilePaths[0])
    else:
        (signalDir, signalFilename) = os.path.split(args.signalsFilePaths)
    med_outputFilePath = os.path.join(signalDir, signalFilename[0:9]+'_med.CPvariant')
    sem_outputFilePath = os.path.join(signalDir, signalFilename[0:9]+'_sem.CPvariant')
    n_outputFilePath = os.path.join(signalDir, signalFilename[0:9]+'_n.CPvariant')

    # Join signals with CPannot
    signals_with_variants = signals.merge(annot,
                                          how='inner',
                                          left_index=True,
                                          right_index=True)

    # Group by and summarize
    signals_by_variants = signals_with_variants.groupby('variant_number')

    signals_by_variants_med = signals_by_variants.median()
    signals_by_variants_sem = signals_by_variants.sem()
    signals_by_variants_n = signals_by_variants.count()[['1']].rename(columns={'1': 'n'})
    
    # Write to output file
    signals_by_variants_med.to_csv(med_outputFilePath, sep='\t')
    signals_by_variants_sem.to_csv(sem_outputFilePath, sep='\t')
    signals_by_variants_n.to_csv(n_outputFilePath, sep='\t')

    return 1

if __name__ == "__main__":
    main()
