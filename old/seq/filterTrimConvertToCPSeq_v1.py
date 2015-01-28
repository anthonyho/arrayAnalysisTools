#!/usr/bin/env python
#
# Script to process barcodes:
# 1. Filter barcodes that doesn't have the correct consensus sequence 
#    immediately after the 16-base barcode
# 2. Trim the barcodes PF into 16 bases
# 3. Output into CPseq format
#
# v1, Anthony Ho, 7/10/2014

## Import libraries
import os, sys
import argparse
import Levenshtein
from numpy import mean


## Phred Q score calculator, default is Phred-33
def qScore(sequence,phred=33):
    if phred == 64:
        return [ord(n)-64 for i,n in enumerate(sequence)]
    else:
        return [ord(n)-33 for i,n in enumerate(sequence)]

## Calculate the average Phred Q score of a given sequence, defualt is Phred-33
def avgQScore(sequence,phred=33):
    if phred == 64:
        l = [ord(n)-64 for i,n in enumerate(sequence)]
    else:
        l = [ord(n)-33 for i,n in enumerate(sequence)]

    return mean(l)


def main():

    ## Get options and arguments from command line
    parser = argparse.ArgumentParser(description="Filter and trim barcodes and output in CPseq format")
    parser.add_argument("length", type=int, help="first n-th bases to be kept")
    parser.add_argument("conSeq", help="consensus sequence")
    parser.add_argument("inputFile", help="input fastq file")
    parser.add_argument("-o", "--outputFile", default="", help="name of the output file")
    parser.add_argument("-q", "--qsCutoff", type=int, default=20, help="cutoff for individual Q score of each base")
    args = parser.parse_args()

    ## Count total number of lines and check if input file is a fastq file
    numLines = 0
    with open(args.inputFile, "r") as r:
        for line in r:
            numLines += 1
            if numLines == 1:
                if line.rstrip().split('\t')[0][0] != "@":
                    print "This script is designed for FASTQ files only!"
                    return 1
    numSeqs = numLines/4
    print "Processing  a total of", str(numSeqs), "sequences..."

    ## Open file for writing
    if args.outputFile == "":
        outputFileName = args.inputFile+".trim.CPseq"
    else:
        outputFileName = args.outputFile+".trim.CPseq"

    ## Initialization
    lineCount = 0 

    lenConSeq = len(args.conSeq)
    endConSeq = lenConSeq + args.length
    
    goodSeqs = 0
    crapSeqs = 0
    lowQSeqs = 0 
    
    ## Going through input file line by line
    with open(args.inputFile, "r") as r, open(outputFileName,"w") as w:
        for line in r:
            # Reading the names, sequence and q scores
            lineCount += 1
            if lineCount % 4 == 1:
                name = line.rstrip().split('\t')[0].split(' ')[0]
            elif lineCount % 4 == 2:
                seq = line.rstrip().split('\t')[0]
            elif lineCount % 4 == 0:
                qs = line.rstrip().split('\t')[0]
                # Business
                if 'N' in seq[0:args.length]: # If there is any N in the barcode
                    crapSeqs += 1
                elif seq[args.length:endConSeq] != args.conSeq: # If the consensus sequence doesn't match
                    crapSeqs += 1
                elif Levenshtein.distance(seq[0:args.length][0]*args.length,seq[0:args.length]) <= 1: # Discard homopolymers 
                    crapSeqs += 1
                elif not(all(i >= args.qsCutoff for i in qScore(qs[0:args.length]))): # If any base has a q score < cutoff
                    lowQSeqs += 1
                else:
                    goodSeqs += 1
                    w.write(name+"\t"+seq[0:args.length]+"\t"+qs[0:args.length]+"\n")
                
    r.close()
    w.close()

    print "Found", str(goodSeqs), "good barcodes. (", str(round(float(goodSeqs)/numSeqs*100,2)), "%)"
    print "Found", str(crapSeqs), "crap barcodes. (", str(round(float(crapSeqs)/numSeqs*100,2)), "%)"
    print "Found", str(lowQSeqs), "low Q-score barcodes. (", str(round(float(lowQSeqs)/numSeqs*100,2)), "%)"

    return 1

if __name__ == "__main__":
    main()
