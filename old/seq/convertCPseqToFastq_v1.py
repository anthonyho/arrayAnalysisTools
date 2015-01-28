#!/usr/bin/env python
#
# Convert CPseq file into fastq files
#
# v1, Anthony Ho, 6/28/2014

## Import libraries
import os, sys
import argparse

def main():

    ## Get options and arguments from command line
    parser = argparse.ArgumentParser(description="Convert CPseq file into fastq files")
    parser.add_argument("PEorSR", choices=["PE","SR"], help="pair-end (PE) or single read (SR) sequences?" )
    parser.add_argument("CPseqFile", help="CPseq file to be converted")
    parser.add_argument("-o", "--outputFile", default="", help="name of the output file")
    args = parser.parse_args()

    ## Determine if the reads are pair-end or single-read 
    if args.PEorSR == "PE":
        PE = True
    else: 
        PE = False
     
    ## Output filename
    if args.outputFile == "":
        r1FilePath = args.CPseqFile+"_r1.fastq"
        if PE:
            r2FilePath = args.CPseqFile+"_r2.fastq"
        i1FilePath = args.CPseqFile+"_i1.fastq"
    else:
        r1FilePath = args.outputFile+"_r1.fastq"
        if PE:
            r2FilePath = args.outputFile+"_r2.fastq"
        i1FilePath = args.outputFile+"_i1.fastq"

    ## Initialization
    lineCount = 0 

    ## Loop through all alignment results in the CPseq file
    if not PE:
        with open(args.CPseqFile, "r") as r, open(r1FilePath, "w") as wr1, open(i1FilePath, "w") as wi1:
            for line in r:
                # Reading line by line
                lineCount += 1
                seqLine = line.rstrip().split('\t')
                name = seqLine[0]
                exp = seqLine[1]
                r1 = seqLine[2]
                Q1 = seqLine[3]
                i1 = seqLine[4]
                Qi = seqLine[5]        
                
                wr1.write("@"+name+"\n")
                wr1.write(r1+"\n")
                wr1.write("+\n")
                wr1.write(Q1+"\n")
            
                wi1.write("@"+name+"\n")
                wi1.write(i1+"\n")
                wi1.write("+\n")
                wi1.write(Qi+"\n")
    else:
        with open(args.CPseqFile, "r") as r, open(r1FilePath, "w") as wr1, open(r2FilePath, "w") as wr2, open(i1FilePath, "w") as wi1:
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
                i1 = seqLine[6]
                Qi = seqLine[7]        
                
                wr1.write("@"+name+"\n")
                wr1.write(r1+"\n")
                wr1.write("+\n")
                wr1.write(Q1+"\n")
                
                wr2.write("@"+name+"\n")
                wr2.write(r2+"\n")
                wr2.write("+\n")
                wr2.write(Q2+"\n")
                
                wi1.write("@"+name+"\n")
                wi1.write(i1+"\n")
                wi1.write("+\n")
                wi1.write(Qi+"\n")

    wr1.close()
    if PE:
        wr2.close()
    wi1.close()

    print "Processed a total of", lineCount, "sequences"

    return 1
 
if __name__ == "__main__":
    main()

