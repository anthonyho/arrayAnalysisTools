#!/usr/bin/env python
#
# Merge pair-end reads together from a CPseq file
#
# v1, Anthony Ho, 6/26/2014

## Import libraries
import os, sys
import argparse
import string
import Levenshtein
#from Bio import SeqIO
#from Bio import AlignIO
import numpy as np

## Reverse complement
complement = string.maketrans("ATCGN", "TAGCN")
def reverseComplement(sequence):
    return sequence.upper().translate(complement)[::-1]

## Align with mismatch, starting from 3' end, find first and move on, assumes only one
def fuzzyAlign(sSeq,lSeq,mismatch,exactOrNot,rangeMax):
    lenSSeq = len(sSeq)
    lenLSeq = len(lSeq)
    # Loop through equal size windows from the 3' end
    # Faster version: scanning for the first alignment result under a given number of mismatches
    if rangeMax != 0:
        mindist = 100
        for i in range(0,lenLSeq-lenSSeq+1)[::-1]:  
            lWindow = lSeq[i:i+lenSSeq]
            dist = Levenshtein.distance(lWindow, sSeq)
            # Find the first match from 3' end then break
            if dist <= rangeMax and dist <= mindist:
                mindist = dist
        if mindist != 100:
            return i, mindist, 1, 0
    # Slower but more accurate version: starting from the lowest number of mismatches first
    if exactOrNot:
        for j in range(rangeMax,mismatch):
            for i in range(0,lenLSeq-lenSSeq+1)[::-1]:  
                lWindow = lSeq[i:i+lenSSeq]
                dist = Levenshtein.distance(lWindow, sSeq)
                # Find the first match from 3' end then break
                if dist == j+1:
                    return i, dist, 0, 1

## Phred Q score calculator, default is Phred-33
def qScore(sequence,phred=33):
    if phred == 64:
        return [ord(n)-64 for i,n in enumerate(sequence)]
    else:
        return [ord(n)-33 for i,n in enumerate(sequence)]

## Calculate the average Phred Q score of a given sequence, defualt is Phred-33
def aveQScore(sequence,phred=33):
    if phred == 64:
        l = [ord(n)-64 for i,n in enumerate(sequence)]
    else:
        l = [ord(n)-33 for i,n in enumerate(sequence)]
    return np.mean(l)

def main():

    ## Get options and arguments from command line
    parser = argparse.ArgumentParser(description="Merge pair-end reads together from a CPseq file")
    parser.add_argument("lenOverlap", type=int, help="length of the overlap used for PE alignment" )
    parser.add_argument("lenEndTrim", type=int, help="the number of bases to be trimmed from the 3' end of r2 before starting the overlap")
    parser.add_argument("mismatchAllowed", type=int, help="maximum number of mismatches within the overlap allowed")
    parser.add_argument("CPseqFile", help="CPseq file containing the PE reads to be merged")
    parser.add_argument("-s", "--speedUp",type=int, default=0, help="")
    parser.add_argument("-o", "--outputFile", default="", help="path of the output file (no suffix)")
    args = parser.parse_args()
     
    ## Output filename
    if args.outputFile == "":
        mergedFilePath = args.CPseqFile+".join"
        unmergedFilePath = args.CPseqFile+".notjoin"
    else:
        mergedFilePath = args.outputFile+".join"
        unmergedFilePath = args.outputFile+".notjoin"

    ## Initialization
    lineCount = 0 
    mismatchStat = [0] * (args.mismatchAllowed + 1)
    fuzzyAlignRangeMax = min(args.speedUp,args.mismatchAllowed)
    exactOrNot = args.mismatchAllowed > args.speedUp
    a=0
    b=0
    c=0

    numCrapSeq = 0
    numLowQSeq = 0


    ## Loop through all alignment results in the CPseq file
    with open(args.CPseqFile, "r") as r, open(mergedFilePath, "w") as wm, open(unmergedFilePath, "w") as wu:
        for line in r:
            # Reading line by line
            lineCount += 1
            seqLine = line.rstrip().split('\t')
            name = seqLine[0]
            exp = seqLine[1]
            r1 = seqLine[2]
            Q1 = seqLine[3]
            r2 = seqLine[4]
            Q2 = seqLine[5]
            bc = seqLine[6]
            Qbc = seqLine[7]        

            if r1[-4:]== "NNNN" or r2[0] == "N":
                numCrapSeq += 1
            elif aveQScore(Q1[-100:]) < 28 or aveQScore(Q2[-100:]) < 28:
                numLowQSeq += 1                             
            else:
                # Align read 2 to read 1
                # Looking for perfect match
                l = len(r2) - args.lenEndTrim 
                r2window = reverseComplement(r2[l-args.lenOverlap:l]) ###!!!
                posMatch = r1.rfind(r2window)
                if posMatch > 0:
                    a+=1
                    mismatchStat[0] += 1
                # If not, allow mismatches
                elif args.mismatchAllowed > 0:
                    alignment = fuzzyAlign(r2window,r1,args.mismatchAllowed,exactOrNot,fuzzyAlignRangeMax)
                    if alignment != None:
                        (posMatch, mismatch, e, f) = alignment
                        mismatchStat[mismatch] += 1
                        b = b+e
                        c = c+f
                        if mismatch == 7:
                            print name
#            wm.write(str(aveQScore(Q1[-100:]))+"\n")
#            wu.write(str(aveQScore(Q2[-100:]))+"\n")
    
    print Qbc
    print qScore (Qbc)
    print aveQScore(Qbc)
    print "Total sequence analyzed:", lineCount
    print "Number of crap sequences:", numCrapSeq
    print "Number of low Q sequences:", numLowQSeq
    print "Mismatch statistics:", mismatchStat

    print a
    print b
    print c

    r.close()
    wm.close()
    wu.close()

    return 1
 
if __name__ == "__main__":
    main()
