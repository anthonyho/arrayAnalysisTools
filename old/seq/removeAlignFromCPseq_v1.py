#!/usr/bin/env python
#
# Remove sequences from a CPseq file with specific 
# SAM flags given the alignment results 
#
# v1, Anthony Ho, 6/26/2014

## Import libraries
import os, sys
import argparse

def main():

    ## Get options and arguments from command line
    parser = argparse.ArgumentParser(description="Remove sequences from a CPseq file with specific SAM flags given the alignment results")
    parser.add_argument("PEorSR", choices=["PE","SR"], help="pair-end (PE) or single read (SR) alignment?" )
    parser.add_argument("alignFile", type=argparse.FileType('r'), help="alignment results in SAM format")
    parser.add_argument("CPseqFile", type=argparse.FileType('r'), help="CPseq file to be modified")
    parser.add_argument("-o", "--outputFile", default="", help="name of the output file")
    flagsGroup = parser.add_mutually_exclusive_group(required=True)
    flagsGroup.add_argument("-k", "--keep", help="alignment flags of the sequences to be kept (multiple flags separated by commas)")
    flagsGroup.add_argument("-r", "--remove", help="alignment flags of the sequences to be removed (multiple flags separated by commas)")
    args = parser.parse_args()

    ## Read SAM flags and determine if user wants to keep or remove sequences with flags
    if args.keep == None:
        k = False
        flags = [x.strip() for x in args.remove.split(',')]
    else: 
        k = True
        flags = [x.strip() for x in args.keep.split(',')]

    ## Determine if the reads are pair-end or single-read 
    if args.PEorSR == "PE":
        pair = 2
    else: 
        pair = 1

    ## Read sequences from the files
    alnList = args.alignFile.readlines()
    seqList = args.CPseqFile.readlines()
 
    ## Determine the number of lines in the SAM header
    i = 0 
    while True:
        if alnList[i].rstrip().split('\t')[0][0] == '@':
            i += 1
        else: 
            break
    numHeaders = i 

    ## Initialization and checking if the number of lines makes sense
    numRemoved = 0
    numAln = len(alnList)
    numSeqs = len(seqList)
    if numSeqs*pair+numHeaders != numAln:
        print "Number of sequences in SAM file does not match with the number of sequences in CPseq file!"
        return 1
    
    ## Open file for writing
    if args.outputFile == "":
        outputFileName = args.CPseqFile.name+".pAln"
    else:
        outputFileName = args.outputFile
    output = open(outputFileName,"w")

    ## Loop through all alignment results in the SAM file
    ## PE:
    if args.PEorSR == "PE":
        for i in range(0,numSeqs):
            alnFlag1 = alnList[i*pair+numHeaders].rstrip().split('\t')[1]
            alnFlag2 = alnList[i*pair+numHeaders+1].rstrip().split('\t')[1]
            if k and (alnFlag1 in flags) and (alnFlag2 in flags):
                output.write(seqList[i])                
            elif not(k) and not(alnFlag1 in flags) and not(alnFlag2 in flags):
                output.write(seqList[i])
            else:
                numRemoved += 1
    ## SR:
    else:
        for i in range(0,numSeqs):
            alnFlag = alnList[i*pair+numHeaders].rstrip().split('\t')[1]
            if k and (alnFlag in flags):
                output.write(seqList[i])                
            elif not(k) and not(alnFlag in flags):
                output.write(seqList[i])
            else:
                numRemoved += 1

    output.close()

    if args.PEorSR =="PE":
        print "\nThe alignment results from the pair-end SAM file:", os.path.abspath(args.alignFile.name) 
    else:
        print "\nThe alignment results from the single-read SAM file:", os.path.abspath(args.alignFile.name) 
    print "is used to remove sequences from the CPseq file:", os.path.abspath(args.CPseqFile.name)
    if k:
        print "without the SAM flags:", flags
    else:
        print "with the SAM flags:", flags
    print "\nTotal number of sequences analyzed:", numSeqs, "(100.0%)"
    print "Number of sequences removed:", numRemoved, "("+str(round(float(numRemoved)/numSeqs*100,2))+"%)"
    print "number of sequences remained:", numSeqs-numRemoved, "("+str(round(float(numSeqs-numRemoved)/numSeqs*100,2))+"%)"

    return 1
 
if __name__ == "__main__":
    main()

