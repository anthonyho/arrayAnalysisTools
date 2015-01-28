#!/usr/bin/env python
#
# Find the corresponding indices (from CPseq) of reads (from fastq) 
# and output in CPseq 
#
# v1, Anthony Ho, 6/28/2014

## Import libraries
import os, sys
import argparse

def main():

    ## Get options and arguments from command line
    parser = argparse.ArgumentParser(description="Find indices of reads")
    parser.add_argument("fastqFile", type=argparse.FileType('r'), help="fastq file with the reads")    
    parser.add_argument("CPseqFile", type=argparse.FileType('r'), help="CPseq file containing the indices")
    parser.add_argument("-o", "--outputFile", default="", help="name of the output file")
    args = parser.parse_args()

    ## Read sequences from the files
    fSeqList = args.fastqFile.readlines()
    cSeqList = args.CPseqFile.readlines()

    ## Output filename
    if args.outputFile == "":
        outputFilePath = args.CPseqFile.name+".CPseq"
    else:
        outputFilePath = args.outputFile.name+".CPseq"
    output = open(outputFilePath,"w")

    ## Initialization
    numFastqLines = len(fSeqList)
    j = 0

    ## Loop through all sequences in the fastq file
    for i in range(0,numFastqLines):
        if i % 4 == 0:
            name = fSeqList[i].rstrip().split('\t')[0][1:]
        elif i % 4 == 1:
            read = fSeqList[i].rstrip().split('\t')[0]
        elif i % 4 == 3:
            qScores = fSeqList[i].rstrip().split('\t')[0]
            foundI = False
            # Find the corresponding index
            while not(foundI):
                cName = cSeqList[j].rstrip().split('\t')[0]
                if cName == name:
                    exp = cSeqList[j].rstrip().split('\t')[1]
                    index = cSeqList[j].rstrip().split('\t')[6]
                    qIndex = cSeqList[j].rstrip().split('\t')[7]
                    foundI = True
                j += 1
            output.write(name+"\t"+exp+"\t"+read+"\t"+qScores+"\t"+index+"\t"+qIndex+"\n")
                
    output.close()

    return 1
 
if __name__ == "__main__":
    main()

