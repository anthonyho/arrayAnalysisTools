#!/usr/bin/env python
#
# Merge and analyze sequences by barcodes
#
# Anthony Ho, ahho@stanford.edu, 4/27/2016
# Last update 8/9/2016


# Import libraries
import os, sys
import argparse
import numpy as np
import pandas as pd
import parlib
import seqlib


# Find the consensus sequence of a barcode block
# Assuming all sequences in a barcode block all have the same length
def consensusVoting(group):

    numSeq = len(group['seq'])
    
    # No need to vote if every sequence agrees
    if group['seq'].tolist() == [group['seq'].iloc[0]] * numSeq:
        consensus = group['seq'].iloc[0]
        numSeqConsent = numSeq
    else:
        consensus = ""
        bases = "ACGTN"
        charArray = np.array(map(list, group['seq']))        
        # Go through all positions
        for i in range(0, len(charArray[0])):
            # Define base array
            baseArray = charArray[:, i].tolist()
            # Count bases and vote
            baseCount = (baseArray.count('A'), 
                         baseArray.count('C'), 
                         baseArray.count('G'), 
                         baseArray.count('T'), 
                         baseArray.count('N'))
            vote = np.argmax(baseCount)
            consensus += bases[vote]
        # Compute barcode block stat
        numSeqConsent = np.sum(group['seq'] == consensus)

    return pd.Series({'seq': consensus,
                      'numSeqConsent': numSeqConsent,
                      'numSeq': numSeq})


def main():

    # Get options and arguments from command line
    parser = argparse.ArgumentParser(description="Merge and analyze sequences by barcodes")
    parser.add_argument('-n', '--numCores', type=int, default=1, help="number of cores to use (default=1)")
    parser.add_argument('-v', '--verbose', type=int, default=0, help="verbosity of progress in parallel mode. 0 = no verbosity. (default=0)")
    parser.add_argument("-q", "--avgQScoreCutoff", type=int, default=28, help="cutoff for average q score of a read")
    parser.add_argument("-b", "--avgQBCScoreCutoff", type=int, default=20, help="cutoff for average q score of a barcode")
    parser.add_argument("inputFilePath", help="path to the CPseq-like file (with header line indicating barcode and seq) to be analyzed")
    parser.add_argument("outputFilePath", help="path to the CPseq-like file (with header line indicating barcode and seq) to be analyzed")
    args = parser.parse_args()

    # Load input file
    print "Loading reads..."
    allReads = pd.read_csv(args.inputFilePath, sep='\t', dtype='object', 
                           usecols=['barcode', 'barcode_Q', 'seq', 'seq_Q'])

    # Filter sequences by their average Q scores
    # skip this step if avgQScoreCutoff and avgQBCScoreCutoff are set to 0
    if args.avgQBCScoreCutoff != 0 or args.avgQScoreCutoff != 0:
        print "Filtering reads by Q scores..."
        if args.numCores == 1:
            bcPF = allReads['barcode_Q'].apply(seqlib.avgQScore) > args.avgQBCScoreCutoff
            readsPF = allReads['seq_Q'].apply(seqlib.avgQScore) > args.avgQScoreCutoff
        else:
            bcPF = parlib.parallelApply(allReads['barcode_Q'], seqlib.avgQScore, args.numCores, args.verbose) > args.avgQBCScoreCutoff
            readsPF = parlib.parallelApply(allReads['seq_Q'], seqlib.avgQScore, args.numCores, args.verbose) > args.avgQScoreCutoff
        filteredReads = allReads[bcPF & readsPF]
    else:
        filteredReads = allReads

    # Group by barcodes and merge in parallel
    print "Grouping and merging sequences with the same barcodes..."
    filteredReads_grouped = filteredReads.groupby('barcode')
    
    if args.numCores == 1:
        mergedReads = filteredReads_grouped.apply(consensusVoting)
        mergedReads.reset_index(inplace=True)
    else:
        mergedReads = parlib.parallelApply(filteredReads_grouped, consensusVoting, args.numCores, args.verbose)
        mergedReads.reset_index(inplace=True)

    # Write to file
    print "Writing to file..."
    mergedReads = mergedReads[['barcode', 'seq', 'numSeqConsent', 'numSeq']]
    mergedReads.to_csv(args.outputFilePath, sep='\t', index=False)

    return 1


if __name__ == "__main__":
    main()
