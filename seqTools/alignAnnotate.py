#!/usr/bin/env python
#
# Align and annotate a given list of sequences
# against a set of reference sequences
#
# Anthony Ho, ahho@stanford.edu, 1/28/2015
# Last update 1/28/2015


## Import libraries
import os, sys
import argparse
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from Bio import SeqIO
import multiprocessing
import seqlib


## Make annotation label given alignment information
def makeAnnotation(label, numMM, numIn, numDel, listMM, listIn, listDel):
    csMM = ','.join(listMM)
    csIn = ','.join(listIn)
    csDel = ','.join(listDel)
    listAnnotation = [label, str(numMM), str(numIn), str(numDel), csMM, csIn, csDel]
    annotation = ':'.join(listAnnotation)
    return annotation


## Align and annotate a query sequence against a set of reference sequences if ignoring indels completely
def alignAnnotateEachSeqMMonly(querySeq, refSeqsDict, startPos, refPosDict, MMcutoff):

    # Make a dictionary of HD of the query against all references that have the same length as query 
    hamDict = {refSeqID: seqlib.hamming(querySeq, refSeq) for refSeqID, refSeq in refSeqsDict.items() if len(querySeq) == len(refSeq)}

    # If any of the references has the same length as query:
    if hamDict:

        # Find the (IDs of) references that have the smallest hamming distance (i.e. most similar to reference)
        minHam = min(hamDict.values())
        minHamRefIDs = [refSeqID for refSeqID, refSeqHam in hamDict.items() if refSeqHam == minHam]
        
        if minHam <= MMcutoff:
            # Output warning if more than one reference has the smallest hamming distance
            if len(minHamRefIDs) > 1:
                print "Warning! Multiple reference sequences have the same Hamming distance with the query sequence"

            # Get list of mismatches
            matchID = minHamRefIDs[0]
            listMM = seqlib.findMismatches(querySeq, refSeqsDict[matchID], startPos, refPosDict[matchID])
            return makeAnnotation(matchID, minHam, 0, 0, listMM, [], [])

        # Report N/A if # MM > MMcutoff
        else:
            return makeAnnotation('NA', np.nan, np.nan, np.nan, [], [], [])

    # Report N/A if indels
    else:
        return makeAnnotation('NA', np.nan, np.nan, np.nan, [], [], [])


## Align and annotate a query sequence against a set of reference sequences considering indels
def alignAnnotateEachSeqMMindels(querySeq, refSeqsDict, startPos, refPosDict):

    # Make a dictionary of alignment results against all references 
    alignmentsDict = {refSeqID: seqlib.alignEmbossNeedle(querySeq, refSeq) for refSeqID, refSeq in refSeqsDict.items()}

    # Make a dictionary of HD of all the aligned query-reference pairs
    hamDict = {refSeqID: seqlib.hamming(alignedQuery, alignedRef) for refSeqID, [alignedQuery, alignedRef] in alignmentsDict.items()}

    # Find the references (ID-seq pairs) that have the smallest hamming distance (i.e. the most similar)
    minHam = min(hamDict.values())
    minHamAlignmentIDs = [refSeqID for refSeqID, refSeqHam in hamDict.items() if refSeqHam == minHam]

    # Output warning if more than one reference has the min hamming distance
    if len(minHamAlignmentIDs) > 1:
        print "Warning! Multiple aligned reference sequences have the same Hamming distance with the aligned query sequence"

    # Get list of mismatches and indels
    matchID = minHamAlignmentIDs[0]
    listMM, listIn, listDel = seqlib.findMismatchesAndIndels(alignmentsDict[matchID][0], alignmentsDict[matchID][1], startPos, refPosDict[matchID])
    return makeAnnotation(matchID, len(listMM), len(listIn), len(listDel), listMM, listIn, listDel)

    
## Align and annotate a query sequence against a set of reference sequences
## Assume all sequences have been converted to upper cases
def alignAnnotateEachSeq(querySeq, refSeqsDict, startPos, refPosDict, MMcutoff, indelMode):
    
    # Check if perfect match with any of the reference sequences
    for refSeqID, refSeq in refSeqsDict.items():
        if querySeq == refSeq:
            return makeAnnotation(refSeqID, 0, 0, 0, [], [], [])

    # Look for mismatches only no-indel mode:
    if not indelMode:
        if not MMcutoff:
            MMcutoff = len(querySeq)
        return alignAnnotateEachSeqMMonly(querySeq, refSeqsDict, startPos, refPosDict, MMcutoff)

    # Look for mismatches and indels in indel mode:
    else:
        return alignAnnotateEachSeqMMindels(querySeq, refSeqsDict, startPos, refPosDict)

def tmpFunc(i, seq, refSeqsDict, startPos, refPosDict, MMcutoff, indelMode):
    print i, alignAnnotateEachSeq(seq, refSeqsDict, startPos, refPosDict, MMcutoff, indelMode)
    return
    

def main():

    ## Get options and arguments from command line
    parser = argparse.ArgumentParser(description="align and annotate a given list of sequences against a set of reference sequences")
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-i', '--indel', action='store_true', help="indel mode: enable alignment and annotation of sequences with indels (default=false)")
    group.add_argument('-m', '--MMcutoff', type=int, help="indicate the maximum number of mismatches allowed in the no-indel mode (default=length of the query sequence)")
    parser.add_argument('-q', '--numberedMode', action='store_true', help="enable non-sequential numbering of the nucleotides along the sequence using the quality score (default=false)")
    parser.add_argument('-s', '--startPos', type=int, default=1, help="number indictating the position of the first base of the read (default=1)")
    parser.add_argument('-n', '--numCore', type=int, default=1, help="number of cores to use (default=1)") 
    parser.add_argument('refSeqFilePath', help="path to the reference sequences file (in FASTA/FASTQ format)")
    parser.add_argument('seqCol', help="column of the file containing the list of sequences to be aligned and annotated (in Python notation)")
    parser.add_argument('seqFilePath', help="path to the file containing the list of sequences to be aligned and annotated")
    parser.add_argument('outCol', help="column of the output file to which the annotation will be written to (in Python notation)")
    parser.add_argument('outputFilePath', help="path to the output file")
    args = parser.parse_args()


    ## Initialization
    phredOffset = 15
    
    ## Read reference sequences (and numbering of nucleotides if in numberedMode)
    refSeqsDict = {}
    refPosDict = {}
    if args.numberedMode:
        for record in SeqIO.parse(args.refSeqFilePath, 'fastq'):
            refSeqsDict[record.id] = str(record.seq.upper())
            refPosDict[record.id] = [i-phredOffset for i in record.letter_annotations['phred_quality']]
    else:
        for record in SeqIO.parse(args.refSeqFilePath, 'fasta'):
            refSeqsDict[record.id] = str(record.seq.upper())
            refPosDict[record.id] = []
        
    ## Load seqfile
    allSeqs = pd.read_csv(args.seqFilePath, sep='\t', header=None)    
    allAnnt = pd.DataFrame(index=allSeqs.index)


    ## Align and annotate
#    if args.numCore == 1:
    allAnnt = allSeqs[int(args.seqCol)].str.upper().apply(alignAnnotateEachSeq, args=(refSeqsDict, args.startPos, refPosDict, args.MMcutoff, args.indel))
    

       
#    print allAnnt

    for i, seq in enumerate(allSeqs[int(args.seqCol)]):
        allAnnt.iat[i,0] = alignAnnotateEachSeq(seq.upper(), refSeqsDict, args.startPos, refPosDict, args.MMcutoff, args.indel)

    
    #> multiprocessing
    #for i, seq in enumerate(allSeqs[int(args.seqCol)]):
#        alignAnnotateEachSeq(seq.upper(), refSeqsDict, args.startPos, refPosDict, args.MMcutoff, args.indel)
    #alignAnnotateEachSeq(seq.upper(), refSeqsDict, args.startPos, refPosDict, args.MMcutoff, args.indel)

#    Parallel(n_jobs=16)(delayed(tmpFunc)(i, seq.upper(), refSeqsDict, args.startPos, refPosDict, args.MMcutoff, args.indel) for i, seq in enumerate(allSeqs[int(args.seqCol)]))

    
    return 1 

if __name__ == "__main__":
    main()
