#!/usr/bin/env python
#
# Trim sequences from a fastq file to the first n-th bases, where n is a parameter given by the user
#
# v1, Anthony Ho, 6/25/2014

## Import libraries
import os, sys
import argparse

def main():

    ## Get options and arguments from command line
    parser = argparse.ArgumentParser(description="Trim sequences to the first n-th bases")
    parser.add_argument("length", type=int, help="first n-th bases to be kept")
    parser.add_argument("inputFile", type=argparse.FileType('r'), help="input fastq file")
    parser.add_argument("-o", "--outputFile", default="", help="name of the output file")
    args = parser.parse_args()

    ## Read sequences from fastq file
    seqList = args.inputFile.readlines()

    ## Check if input file is a fastq file
    if seqList[0][0] != "@":
        print "This script is designed for FASTQ files only!"
        return 1

    ## Open file for writing
    if args.outputFile == "":
        outputFileName = args.inputFile.name+".trim"
    else:
        outputFileName = args.outputFile
    output = open(outputFileName,"w")

    ## Initialization
    numSeqs = len(seqList)

    ## Loop through all sequences
    ## and trim lengths of the sequences
    for i in range(0,numSeqs):
        if i % 2 == 0:
            output.write(seqList[i])           
        else:
            output.write(seqList[i][0:args.length]+'\n')

    output.close()

    return 1

if __name__ == "__main__":
    main()
