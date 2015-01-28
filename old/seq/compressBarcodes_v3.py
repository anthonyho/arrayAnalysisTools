#!/usr/bin/env python
#
# Group and analyze sequences by barcode blocks
#
# Anthony Ho, 6/29/2014
# Based on Jason Buenrostro's script pyCompressBBsv2.py

## Import libraries
import os, sys
import argparse
import numpy as np
import Levenshtein
from collections import Counter

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
    return np.mean(l)


## Finding the consensus sequence of a barcode block
## !!! Needs to be refined !!!
def consensusVoting(r1Block,Q1Block,degeneracy):
    ## Find the consensus sequence
    consensus = ""
    charArray = np.array(map(list, r1Block[:]))
    bases = "ACGT"
    for i in range(0,len(charArray[0])):
        ## Base array
        baseArray = charArray[:,i].tolist()
        ## Count bases and vote
        baseCount = (baseArray.count('A'), baseArray.count('C'), baseArray.count('G'), baseArray.count('T'))
        vote = np.argmax(baseCount)
        consensus += bases[vote]
    return consensus


def main():

    ## Get options and arguments from command line
    parser = argparse.ArgumentParser(description="Group and analyze sequences by barcode blocks")
    parser.add_argument("CPseqFile", help="CPseq file to be analyzed")
    parser.add_argument("-q", "--avgQScoreCutoff", type=int, default=28, help="cutoff for average q score of the read")
    parser.add_argument("-b", "--avgQBCScoreCutoff", type=int, default=20, help="cutoff for average q score of the barcode")
    parser.add_argument("-d", "--degeneracyCutoff", type=int, default=2, help="cutoff for the degeneracy of a barcode block")
    ## Should be changed to 3
#    parser.add_argument("-c", "--conSeq1", default="", help="constant sequence")
#    parser.add_argument("-n", "--corrN", action="store_true", help="corrects N bases")
    args = parser.parse_args()

    ## Open file for writing
    outputFileName = args.CPseqFile+".unq"
    statFileName = args.CPseqFile+".stat"

    ## Initialization
    numSeqs = 0 
    lineCount = 0

    r1Block = []
    Q1Block = []
    aQbBlock = []

    # Stat on reads
    goodSeqs = 0
    crapSeqs = 0
    uncertainBCSeqs = 0
    badBCSeqs = 0
    
    # Stat of barcodes
    goodBC = 0
    uncertainBC = 0
    badBC = 0

    # Stat for merging
    perfectSeqs = 0
    sameLenSeqs = 0
    diffLenSeqs = 0

    ## Initial counting
    with open(args.CPseqFile,"r") as r:
        for line in r:
            numSeqs += 1
            if numSeqs == 1:
                lastBC = line.rstrip().split('\t')[4]

    ## Going through the CPseq file
    with open(args.CPseqFile, "r") as r, open(outputFileName, "w") as wo, open(statFileName, "w") as ws:
        for line in r:

            # Reading line by line
            lineCount += 1
            if lineCount % 10000 == 0:
                print "Processing the "+str(lineCount) +"th sequence" 
            seqLine = line.rstrip().split('\t')
            r1 = seqLine[2]
            Q1 = seqLine[3]
            BC = seqLine[4]
            Qbc = seqLine[5]

            # At the beginning of each barcode block, 
            # and do it for the very last line instead of the very first line
            if BC != lastBC or lineCount == numSeqs:
                
                # Add in the very last line to the barcode block
                # Append sequence and q-scores to the barcode block
                # Check the sequences before putting into block
                if lineCount == numSeqs:
                    if avgQScore(Q1) <= args.avgQScoreCutoff:
                        crapSeqs += 1
                    else:
                        r1Block.append(r1)
                        Q1Block.append(Q1)
                        aQbBlock.append(avgQScore(Qbc))

                # Analyze barcode block here
                degeneracy = len(r1Block)
                if degeneracy > 0:
                    r1Block = np.array(r1Block)
                    Q1Block = np.array(Q1Block)
                    aQbBlock = np.array(aQbBlock)
                    
                    # Discard homopolymer or near-homopolymer barcodes
                    if np.mean(aQbBlock) < args.avgQBCScoreCutoff:
                        badBC += 1
                        badBCSeqs += degeneracy
                    elif Levenshtein.distance(lastBC[0]*len(lastBC),lastBC) <= 1:
                        badBC += 1
                        badBCSeqs += degeneracy
                    # Discard barcode block with degeneracy <= 3
                    elif degeneracy <= args.degeneracyCutoff: 
                        uncertainBC += 1
                        uncertainBCSeqs += degeneracy
                    # Otherwise proceed to merging and statistics
                    else: 
                        # Dealing with statistics
                        ws.write(str(degeneracy)+"\n")
                        goodBC += 1
                        goodSeqs += degeneracy
                        
                        # Merging
                        # !!! Needs to be refined !!!
                        # Find perfect matches within a block first:
                        if r1Block[:].tolist() == [r1Block[0]]*degeneracy:
                            perfectSeqs += 1
                            wo.write(r1Block[0]+"\t"+lastBC+"\t"+str(degeneracy)+"\t0\t"+str(len(r1Block[0]))+"\n")
                        # Then deal with sequence blocks of the same length of sequences:
                        elif max(r1Block,key=len) == min(r1Block,key=len):
                            sameLenSeqs += 1
                            consensusSeq = consensusVoting(r1Block,Q1Block,degeneracy)
                            wo.write(consensusSeq+"\t"+lastBC+"\t"+str(degeneracy)+"\t1\t"+str(len(r1Block[0]))+"\n")
                        # If not, remove the sequences that are not of the same length:
                        else:
                            lenList = [len(x) for x in r1Block]
                            lenCount = Counter(lenList)
                            modeLen = lenCount.most_common(1)[0][0]
                            truncatedR1Block = [x for x in r1Block if len(x) == modeLen]
                            diffLenSeqs += 1
                            consensusSeq = consensusVoting(truncatedR1Block,Q1Block,degeneracy)
                            wo.write(consensusSeq+"\t"+lastBC+"\t"+str(degeneracy)+"\t2\t"+str(modeLen)+"\n")
                
                # Initialize for the next barcode block
                r1Block = []
                Q1Block = []
                aQbBlock = []

            # Append sequence and q-scores to the barcode block
            # Check the sequences before putting into block
            if avgQScore(Q1) <= args.avgQScoreCutoff:
                crapSeqs += 1
            else:
                r1Block.append(r1)
                Q1Block.append(Q1)
                aQbBlock.append(avgQScore(Qbc))

            # Make the current barcode the new last barcode
            lastBC = BC
               
    ## Printing summary
    print "\nTotal number of sequences analyzed:", str(numSeqs), "(100%)"
    print "Number of sequences passing filter:", str(goodSeqs), "("+str(round(float(goodSeqs)/numSeqs*100,2)), "%)"
    print "Number of low quality sequences:", str(crapSeqs), "("+str(round(float(crapSeqs)/numSeqs*100,2)), "%)"
    print "Number of sequences with bad BC:", str(badBCSeqs), "("+str(round(float(badBCSeqs)/numSeqs*100,2)), "%)"
    print "Number of sequences with not enough BC degeneracy:", str(uncertainBCSeqs), "("+str(round(float(uncertainBCSeqs)/numSeqs*100,2)), "%)"

    totalNumBC = goodBC + uncertainBC + badBC
    print "\nTotal number of barcodes passing filter:", str(totalNumBC), "(100%)"
    print "Number of unique sequences:", str(goodBC), "("+str(round(float(goodBC)/totalNumBC*100,2)), "%)"
    print "Number of bad barcodes:", str(badBC), "("+str(round(float(badBC)/totalNumBC*100,2)), "%)"
    print "Number of barcodes with not enough BC degeneracy:", str(uncertainBC), "("+str(round(float(uncertainBC)/totalNumBC*100,2)), "%)"

    ## Testing
    print "\nNumber of unique sequence with same sequences within a barcode block:", str(perfectSeqs), "("+str(round(float(perfectSeqs)/goodBC*100,2)), "%)"
    print "Number of unique sequences with different sequences of all the same length within a barcode block:", str(sameLenSeqs), "("+str(round(float(sameLenSeqs)/goodBC*100,2)), "%)"
    print "Number of unique sequences with different sequences of different length", str(diffLenSeqs), "("+str(round(float(diffLenSeqs)/goodBC*100,2)), "%)"

    r.close()
    wo.close()
    ws.close()
    
    return 1

if __name__ == "__main__":
    main()

