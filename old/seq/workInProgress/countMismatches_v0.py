#!/usr/bin/env python
#
# Count mismatches from count*txt files
#
# v0, Anthony Ho, 6/30/2014

## Import libraries
import os, sys
import argparse

def main():

    ## Get options and arguments from command line
    parser = argparse.ArgumentParser(description="Count mismatch from count*txt file")
    parser.add_argument("refFastq", type=argparse.FileType('r'), help="reference fasta file")
    parser.add_argument("seqFile", help="sequences to be counted")
    parser.add_argument("-o", "--outputFile", default="", help="name of the output file")
    args = parser.parse_args()

    ## Read reference fastq file
    refList = args.refFastq.readlines()
    ref = refList[1]

    ## Output filename
    if args.outputFile == "":
        outputFilePath = args.seqFile+".stat"
    else:
        outputFilePath = args.seqFile+".stat"
    output = open(outputFilePath,"w")

    ## Initialization
    lineCount = 0 
    sameLenSeqs = 0
    wrongLenSeqs = 0
    MMList = []
    
    ## Loop through all sequences and count mismatches
    with open(args.seqFile, "r") as r:
        for line in r:
            # Reading line by line
            lineCount += 1
            seqLine = line.rstrip().split('\t')
            seq = seqLine[0]
            bc = seqLine[1]

            if len(seq) != len(ref):
                wrongLenSeqs += 1
            else:
                sameLenSeqs += 1
                MM = 0
                for i in range(0,len(ref)):
                    if seq[i] != ref[i]:
                        MM += 1
                MMList.append(MM)

    for i in range(0,len(MMList)):
        output.write(str(MMList[i])+"\n")

    r.close()
    output.close()

    print "Processed a total of", lineCount, "sequences"
    print "Counted", sameLenSeqs, "sequences with the right length"
    print "Counted", wrongLenSeqs, "sequences with the wrong length"

    return 1
 
if __name__ == "__main__":
    main()

