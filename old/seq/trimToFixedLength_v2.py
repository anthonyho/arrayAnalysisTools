#!/usr/bin/env python
#
# Trim sequences from a fastq file to the first n-th bases, where n is a parameter given by the user
#
# v2, Anthony Ho, 9/30/2014

## Import libraries
import os, sys
import argparse

def main():

    ## Get options and arguments from command line
    parser = argparse.ArgumentParser(description="Trim sequences to the first n-th bases")
    parser.add_argument("length", type=int, help="the first n-th bases to be kept")
    parser.add_argument("inputFile", help="path to the input FASTQ file")
    parser.add_argument("-o", "--outputFile", default="", help="path to the output file")
    args = parser.parse_args()

    ## Path to the output file
    if args.outputFile == "":
        outputFileName = args.inputFile+".trim"
    else:
        outputFileName = args.outputFile

    ## Check if input file is a fastq file before proceeding
    with open(args.inputFile, 'r') as r:
        firstLine = r.readline()
        if firstLine[0] != "@":
            print "This script is designed for FASTQ files only!"
            return 1

    ## Initialization
    lineCount = 0 

    ## Going through the FASTQ file line by line 
    ## to trim the sequences and Q scores
    with open(args.inputFile, "r") as r, open(outputFileName, "w") as w:
        for line in r:
            
            ## Update line count
            lineCount += 1
            
            ## Determine if this line is a head/sequence/spacer/Q score
            ## and write to the output file appropriately
            if lineCount % 2 != 0:
                w.write(line)
            else:
                w.write(line[0:args.length]+'\n')

    ## Printing summary
    print "Total number of sequences processed:", str(lineCount/4)
    print "Trimmed sequences written to:", outputFileName

    r.close()
    w.close()

    return 1

if __name__ == "__main__":
    main()
