#!/usr/bin/env python
#
# Count number of mutations from reads
#
# v1, Anthony Ho, 7/9/2014

## Import libraries
import os, sys
import argparse

## Hamming distance calculator, assuming equal length
def hammingDistance(s1,s2):
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

## Stupid way to find the first (and only) deletion
def findIndel(q,refU):
    for i in range(0,len(refU)):
        if q[i] != refU[i]:
            bail = hammingDistance(q[i+1:i+9],refU[i+2:i+10]) <= 2
            if bail:
                return i+1
    return -1000

def main():

    ## Get options and arguments from command line
    parser = argparse.ArgumentParser(description="Count numbers of mutations from reads")
    parser.add_argument("findIndels", choices=['1','0'], help="find the positions of the indels or not")
    parser.add_argument("mmCutoff", type=int, help="Up to and below which the number of mutation the reads will be output into a file")
    parser.add_argument("refFastq", type=argparse.FileType('r'), help="reference fastq file")
    parser.add_argument("seqFile", help="sequences to be counted")
    parser.add_argument("-o", "--outputFile", default="", help="name of the output file")
    args = parser.parse_args()

    ## Read reference fastq file
    refList = args.refFastq.readlines()
    # Check if fastq file
    if refList[0][0] != '@':
        print "Input file must be fastq file!!"
        return 1
    # Check number of lines in fastq file
    elif len(refList) != 2:
        print "Input file must be fastq file with one reference sequence!!"
        return 1
    # Identify all mutation positions and produce strings of mutation sequence and consensus sequence
    else:
        ref = refList[1] # Get reference sequence
        mutPos = [i for i,N in enumerate(ref) if ref[i] != ref.upper()[i]] # Find the position of the bases that can be mutated
        refCon = ''.join([N for i,N in enumerate(ref) if not(i in mutPos)]) # Make a string of the consensus bases
        numMutPos = len(mutPos)
        refU = ref.upper()

    ## Open files for writing
    mmOutputFilePaths = []
    mmOutput = []
    if args.outputFile == "":
        statOutputFilePath = args.seqFile+".stat"
        if args.findIndels:
            indelsOutputFilePath = args.seqFile+".indel"
        for i in range(0,args.mmCutoff+1):
            mmOutputFilePaths.append(args.seqFile+".mm"+str(i))
    else:
        statOutputFilePath = args.outputFile+".stat"
        if args.findIndels:
            indelsOutputFilePath = args.outputFile+".indel"
        for i in range(0,args.mmCutoff+1):
            mmOutputFilePaths.append(args.outputFile+".mm"+str(i))
    statOutput = open(statOutputFilePath,"w")
    if args.findIndels:
        indelsOutput = open(indelsOutputFilePath,"w")
    for i in range(0,args.mmCutoff+1):
        mmOutput.append(open(mmOutputFilePaths[i],"w"))

    ## Initialization
    refLen = len(refU)

    lineCount = 0 

    sameLenSeqs = 0
    wrongLenSeqs = 0
    badBCSeqs = 0

    correctMutSeqs = 0
    wrongMutSeqs = 0

    MMList = []
    indelList = []
    
    ## Loop through all sequences and count mismatches
    with open(args.seqFile, "r") as r:
        for line in r:
            # Reading line by line
            lineCount += 1
            seqLine = line.rstrip().split('\t')
            seq = seqLine[0]
            bc = seqLine[1]

            # Discard sequences with barcode containing any N
            if 'N' in bc:
                badBCSeqs += 1
            # If sequence is not of the right length
            elif len(seq) != refLen:
                wrongLenSeqs += 1
                # Analyze the n-1 fragments
                if args.findIndels and (len(seq)+1) == refLen:
                    delPos = findIndel(seq,refU)
                    if delPos > 0:
                        indelList.append(delPos)
            # Business
            else:
                sameLenSeqs += 1
                seqCon = ''.join([N for i,N in enumerate(seq) if not(i in mutPos)])
                # If everything outside of the mutation zones matches 
                if seqCon == refCon:
                    correctMutSeqs += 1
                    MM = 0 
                    for i in mutPos:
                        if seq[i] != refU[i]:
                            MM += 1
                    MMList.append(MM)
                    if MM <= args.mmCutoff:
                        seqMut = ''.join([N for i,N in enumerate(seq) if i in mutPos])
                        mmOutput[MM].write(seqMut+"\t"+bc+"\n")
                else:
                    wrongMutSeqs += 1

    ## Writing stat file
    for i in range(0,len(MMList)):
        statOutput.write(str(MMList[i])+"\n")
    ## Writing indel file
    for i in range(0,len(indelList)):
        indelsOutput.write(str(indelList[i])+"\n")

    r.close()
    statOutput.close()
    indelsOutput.close()

    print "\nProcessed a total of", lineCount, "sequences"
    print "Counted", badBCSeqs, "sequences with the badBC"
    print "Counted", wrongLenSeqs, "sequences with the wrong length"
    print "Counted", sameLenSeqs, "sequences with the right length"
    print "\nOf which, "
    print wrongMutSeqs, "are sequences with mutation in the wrong place"
    print correctMutSeqs, "are sequences with mutations in the right place"

    return 1
 
if __name__ == "__main__":
    main()

