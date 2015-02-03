#!/usr/bin/env python
#
# Alignemnt and annotation tools
#
# Anthony Ho, ahho@stanford.edu, 1/28/2015
# Last update 1/31/2015


## Import libraries
import subprocess
import numpy as np

## Parse fasta sequences from string, separated by \n
def parseFastaSeqsAsString(fastaStr, IDchar='>'):
    return [''.join(record.split('\n')[1:]) for record in fastaStr.split(IDchar)[1:]]

## Compute Hamming distance between two sequences, assuming equal length
def hamming(s1, s2):
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

## Find mismatches between two sequences, assuming equal length
def findMismatches(query, reference, startPos=1, refPos=[]):
    # User-defined nucleotide positions
    if not refPos:
        realRefPos = range(0,len(query))
    else:
        realRefPos = np.cumsum(refPos)-1
    # Return list of mismatches
    return [chR+str(startPos+realRefPos[i])+chQ for i, (chQ, chR) in enumerate(zip(query, reference)) if chQ != chR]

## Find mismatches and indels between two aligned sequences (i.e. equal length after gapping)
def findMismatchesAndIndels(query, reference, startPos=1, refPos=[]):

    # User-defined nucleotide positions
    if not refPos:
        realRefPos = range(0,len(query+1))
    else:
        realRefPos = np.cumsum(refPos)-1

    # Find insertions
    # Get all the insertions in the unadjusted coordinates
    rawTupleListIn = [(i, chQ) for i, (chQ, chR) in enumerate(zip(query, reference)) if chR == '-']
    # Adjust to the correct coordinates, map to the user-specific position, and output list of annotations
    listIn = ['+'+str(startPos+realRefPos[i-j])+chQ for j, (i, chQ) in enumerate(rawTupleListIn)]
    
    # Shorten the sequences to get rid of the gaps in the reference sequences and the corresponding 
    # nucleotide in the query sequence
    rawListInPos = [i for i, chQ in rawTupleListIn]
    shortQuery = ''.join([ch for i, ch in enumerate(query) if i not in rawListInPos])
    shortReference = ''.join([ch for i, ch in enumerate(reference) if i not in rawListInPos])
    
    # Find mismatches and deletions
    listMM = [chR+str(startPos+realRefPos[i])+chQ for i, (chQ, chR) in enumerate(zip(shortQuery, shortReference)) if chQ != chR and chQ != '-']
    listDel = [chR+str(startPos+realRefPos[i])+'x' for i, (chQ, chR) in enumerate(zip(shortQuery, shortReference)) if chQ == '-']

    return listMM, listIn, listDel

## Align a query sequence to a reference sequence using EMBOSS's Needleman-Wunsch program
## Note: it assumes the EMBOSS needle program is installed and only works on Unix/Mac-based OS 
def alignEmbossNeedle(query, reference, gapopen=10, gapextend=0.5):
    cmdBash = ' /bin/bash -c '
    cmdAlign = ' needle -asequence <(echo ' + query + ') -bsequence <(echo ' + reference + ')' 
    optionsAlign = ' -gapopen ' + str(gapopen) + ' -gapextend ' + str(gapextend) + ' -aformat fasta -stdout -auto'
    fullCmd = cmdBash + '\"' + cmdAlign + optionsAlign + '\"'
    process = subprocess.Popen(fullCmd, stdin=subprocess.PIPE, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
    process.wait()
    [stdout, stderr] = process.communicate()
    if stderr:
        print 'Needleman-Wunsch alignement returns with error'
    return parseFastaSeqsAsString(stdout)
