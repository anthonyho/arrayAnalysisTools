# Anthony Ho, ahho@stanford.edu, 1/28/2015
# Last update 3/10/2015
"""Library of sequence analysis, alignemnt and annotation tools"""


import subprocess
import numpy as np
import pandas as pd
import string
from itertools import product


# Get a list of annotations of saturated mutations in all given positions
def getSatMutations(listSeqPos, consensus='', rna=False):

    # Make default list of all bases
    if rna:
        allBases = ['A', 'C', 'G', 'U']
    else:
        allBases = ['A', 'C', 'G', 'T']
    # Extract the list for positions from listSeqPos for sorting
    listPos = [seqPos[1:] for seqPos in listSeqPos]
    # Sort listPos and listSeqPos by listPos
    listPos, listSeqPos = zip(*sorted(zip(listPos, listSeqPos)))
    # Get a list of list of all single mutations
    listListSingleMutations = [[seqPos+base for base in allBases] for seqPos in listSeqPos]
    # Get a list of tuple of all possible mutations
    listTupleSatMutations = list(product(*listListSingleMutations))

    # Convert the list of tuple of all possible mutations into
    # a list of annotations of all possible mutations
    listAnnotations = []
    for var in listTupleSatMutations:
        mutations = [mut for mut in var if mut[0] != mut[-1]]
        numMut = len(mutations)
        annotation = consensus + ':' + str(numMut) + ':0:0:' + ','.join(mutations) + ':::'
        listAnnotations.append(annotation)

    return pd.DataFrame({'mutations': listTupleSatMutations,
                         'annotation': listAnnotations})


# Get all four bases as a list
def allBases(rna=False):
    if rna:
        return ['A', 'C', 'G', 'U']
    else:
        return ['A', 'C', 'G', 'T']


# Get all other bases as a list
def allOtherBases(bases, rna=False):
    if rna:
        allBases = ['A', 'C', 'G', 'U']
    else:
        allBases = ['A', 'C', 'G', 'T']
    for base in bases:
        base = base.upper()
        if base in allBases:
            allBases.remove(base)
    return allBases


# Find reverse complement
def reverseComplement(seq, rna=False):
    if rna:
        complements = string.maketrans("AUCGN", "UAGCN")
        return convertToRNA(seq).translate(complements)[::-1]
    else:
        complements = string.maketrans("ATCGN", "TAGCN")
        return convertToDNA(seq).translate(complements)[::-1]


# Standardize the sequence to all upper case and to the correct DNA/RNA base
def standardize(seq, rna=False):
    if rna:
        return convertToRNA(seq)
    else:
        return convertToDNA(seq)
    

# Convert RNA sequences to DNA
def convertToDNA(seq):
    return seq.upper().replace('U', 'T')


# Convert DNA sequences to RNA
def convertToRNA(seq):
    return seq.upper().replace('T', 'U')


# Compute Phred Q score, default is Phred-33
def qScore(seq, phred=33):
    if phred == 64:
        return [ord(n)-64 for i, n in enumerate(seq)]
    else:
        return [ord(n)-33 for i, n in enumerate(seq)]


# Compute the average Phred Q score of a given sequence, defualt is Phred-33
def avgQScore(seq, phred=33):
    if phred == 64:
        l = [ord(n)-64 for i, n in enumerate(seq)]
    else:
        l = [ord(n)-33 for i, n in enumerate(seq)]
    return np.mean(l)


# Parse fasta sequences from string, separated by \n
def parseFastaSeqsAsString(fastaStr, IDchar='>'):
    return [''.join(record.split('\n')[1:]) for record in fastaStr.split(IDchar)[1:]]


# Compute Hamming distance between two sequences, assuming equal length
def hamming(s1, s2):
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


# Find mismatches between two sequences, assuming equal length
def findMismatches(query, reference, startPos=1, refPos=[]):

    # User-defined nucleotide positions
    if not refPos:
        realRefPos = range(0, len(query))
    else:
        realRefPos = np.cumsum(refPos)-1

    # Return list of mismatches
    listMM = [chR+str(startPos+realRefPos[i])+chQ for i, (chQ, chR) in enumerate(zip(query, reference)) if chQ != chR]
    return listMM


# Find mismatches and indels between two aligned sequences (i.e. equal length after gapping)
def findMismatchesAndIndels(query, reference, startPos=1, refPos=[]):

    # User-defined nucleotide positions
    if not refPos:
        realRefPos = range(0, len(query+1))
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


# Align a query sequence to a reference sequence using EMBOSS's Needleman-Wunsch program
# Note: it assumes the EMBOSS needle program is installed and only works on Unix/Mac-based
# OS due to the bash interface
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
