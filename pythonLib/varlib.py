# Anthony Ho, ahho@stanford.edu, 10/1/2015
# Last update 10/1/2015
"""Library functions to generate various sequence variants in library design"""


#import subprocess
#import numpy as np
#import pandas as pd
import seqlib
#import string
#from itertools import product


# Generate a list of variants with single bulges
def generate_1x0(consensus, rna=False):

    seq = seqlib.standardize(consensus, rna)
    listVariants = []
    for i, base in enumerate(seq[:-1]):
        for newBase in seqlib.allBases(rna):
            variant = seq[:i+1] + newBase + seq[i+1:]
            listVariants.append(variant)

    return listVariants

# Generate a list of variants with single mismatches
def generate_1x1(consensus, rna=False):
    
    seq = seqlib.standardize(consensus, rna)
    listVariants = []
    for i, base in enumerate(seq):
        for otherBase in seqlib.allOtherBases(base, rna):
            variant = seq[:i] + otherBase + seq[i+1:]
            listVariants.append(variant)
    
    return listVariants

# Generate a list of variants with single deletions
def generate_0x1(consensus, rna=False):
    
    seq = seqlib.standardize(consensus, rna)
    listVariants = []
    for i, base in enumerate(seq):
        variant = seq[:i] + seq[i+1:]
        listVariants.append(variant)
    
    return listVariants

###########################################################

# Generate a list of variants with 2x0 bulges
def generate_2x0(consensus, rna=False):
    
    seq = seqlib.standardize(consensus, rna)
    listVariants = []
    for i, base in enumerate(seq[:-1]):
        for newBase1 in seqlib.allBases(rna):
            for newBase2 in seqlib.allBases(rna):
                variant = seq[:i+1] + newBase1 + newBase2 + seq[i+1:]
                listVariants.append(variant)
    
    return listVariants

# Generate a list of variants with 2x1 internal loops
def generate_2x1(consensus, rna=False):
    
    seq = seqlib.standardize(consensus, rna)
    listVariants = []
    for i, base in enumerate(seq):
        for otherBase1 in seqlib.allOtherBases(base, rna):
            for otherBase2 in seqlib.allOtherBases(base, rna):
                variant = seq[:i] + otherBase1 + otherBase2 + seq[i+1:]
                listVariants.append(variant)
    
    return listVariants

# Generate a list of variants with 2x2 internal loops
def generate_2x2(consensus, rna=False):
    
    seq = seqlib.standardize(consensus, rna)
    listVariants = []
    for i, base in enumerate(seq[:-1]):
        nextBase = seq[i+1]
        for otherBase1 in seqlib.allOtherBases(base, rna):
            for otherBase2 in seqlib.allOtherBases(nextBase, rna):
                variant = seq[:i] + otherBase1 + otherBase2 + seq[i+2:]
                listVariants.append(variant)
    
    return listVariants

# Generate a list of variants with 1x2 internal loops
def generate_1x2(consensus, rna=False):
    
    seq = seqlib.standardize(consensus, rna)
    listVariants = []
    for i, base in enumerate(seq[:-1]):
        nextBase = seq[i+1]
        for otherBase in seqlib.allOtherBases([base, nextBase], rna):
            variant = seq[:i] + otherBase + seq[i+2:]
            listVariants.append(variant)
    
    return listVariants

# Generate a list of variants with 0x2 deletions
def generate_0x2(consensus, rna=False):
    
    seq = seqlib.standardize(consensus, rna)
    listVariants = []
    for i, base in enumerate(seq[:-1]):
        variant = seq[:i] + seq[i+2:]
        listVariants.append(variant)
    
    return listVariants

###########################################################

# Generate a list of variants with 3x0 bulges
def generate_3x0(consensus, rna=False):
    
    seq = seqlib.standardize(consensus, rna)
    listVariants = []
    for i, base in enumerate(seq[:-1]):
        for newBase1 in seqlib.allBases(rna):
            for newBase2 in seqlib.allBases(rna):
                for newBase3 in seqlib.allBases(rna):
                    variant = seq[:i+1] + newBase1 + newBase2 + newBase3 + seq[i+1:]
                    listVariants.append(variant)
    
    return listVariants

# Generate a list of variants with 3x1 internal loops
def generate_3x1(consensus, rna=False):
    
    seq = seqlib.standardize(consensus, rna)
    listVariants = []
    for i, base in enumerate(seq):
        for newBase1 in seqlib.allBases(rna):
            for newBase2 in seqlib.allBases(rna):
                for newBase3 in seqlib.allBases(rna):
                    variant = seq[:i] + newBase1 + newBase2 + newBase3 + seq[i+1:]
                    listVariants.append(variant)
    
    return listVariants

# Generate a list of variants with 3x2 internal loops
def generate_3x2(consensus, rna=False):
    
    seq = seqlib.standardize(consensus, rna)
    listVariants = []
    for i, base in enumerate(seq[:-1]):
        for newBase1 in seqlib.allBases(rna):
            for newBase2 in seqlib.allBases(rna):
                for newBase3 in seqlib.allBases(rna):
                    variant = seq[:i] + newBase1 + newBase2 + newBase3 + seq[i+2:]
                    listVariants.append(variant)
    
    return listVariants

# Generate a list of variants with 3x3 internal loops
def generate_3x3(consensus, rna=False):
    
    seq = seqlib.standardize(consensus, rna)
    listVariants = []
    for i, base in enumerate(seq[:-2]):
        nextBase = seq[i+1]
        nextNextBase = seq[i+2]
        for otherBase1 in seqlib.allOtherBases(base, rna):
            for newBase2 in seqlib.allBases(rna):
                for otherBase3 in seqlib.allOtherBases(nextNextBase, rna):
                    variant = seq[:i] + otherBase1 + newBase2 + otherBase3 + seq[i+3:]
                    listVariants.append(variant)
    
    return listVariants

# Generate a list of variants with 2x3 internal loops
def generate_2x3(consensus, rna=False):
    
    seq = seqlib.standardize(consensus, rna)
    listVariants = []
    for i, base in enumerate(seq[:-2]):
        nextBase = seq[i+1]
        nextNextBase = seq[i+2]
        for otherBase1 in seqlib.allOtherBases(base, rna):
            for otherBase3 in seqlib.allOtherBases(nextNextBase, rna):
                variant = seq[:i] + otherBase1 + otherBase3 + seq[i+3:]
                listVariants.append(variant)
    
    return listVariants

# Generate a list of variants with 1x3 internal loops
def generate_1x3(consensus, rna=False):
    
    seq = seqlib.standardize(consensus, rna)
    listVariants = []
    for i, base in enumerate(seq[:-2]):
        for newBase in seqlib.allBases(rna):
            variant = seq[:i] + newBase + seq[i+3:]
            listVariants.append(variant)
    
    return listVariants

# Generate a list of variants with 0x3 internal loops
def generate_0x3(consensus, rna=False):
    
    seq = seqlib.standardize(consensus, rna)
    listVariants = []
    for i, base in enumerate(seq[:-2]):
        variant = seq[:i] + seq[i+3:]
        listVariants.append(variant)
    
    return listVariants

###########################################################

# Generate a list of variants with 4x0 bulges
def generate_4x0(consensus, rna=False):
    
    seq = seqlib.standardize(consensus, rna)
    listVariants = []
    for i, base in enumerate(seq[:-1]):
        for newBase1 in seqlib.allBases(rna):
            for newBase2 in seqlib.allBases(rna):
                for newBase3 in seqlib.allBases(rna):
                    for newBase4 in seqlib.allBases(rna):
                        variant = seq[:i+1] + newBase1 + newBase2 + newBase3 + newBase4 + seq[i+1:]
                        listVariants.append(variant)
    
    return listVariants

# Generate a list of variants with 4x4 internal loops
def generate_4x4(consensus, rna=False):
    
    seq = seqlib.standardize(consensus, rna)
    listVariants = []
    for i, base in enumerate(seq[:-3]):
        nextBase = seq[i+1]
        nextNextBase = seq[i+2]
        nextNextNextBase = seq[i+3]
        for otherBase1 in seqlib.allOtherBases(base, rna):
            for newBase2 in seqlib.allBases(rna):
                for newBase3 in seqlib.allBases(rna):
                    for otherBase4 in seqlib.allOtherBases(nextNextNextBase, rna):
                        variant = seq[:i] + otherBase1 + newBase2 + newBase3 + otherBase4 + seq[i+4:]
                        listVariants.append(variant)
    
    return listVariants

###########################################################

# Generate a list of variants with double mismatches
def generate_1x1_1x1(consensus, rna=False):
    
    seq = seqlib.standardize(consensus, rna)
    listVariants = []
    for i, base1 in enumerate(seq):
        for j, base2 in enumerate(seq[i+1:]):
            for otherBase1 in seqlib.allOtherBases(base1, rna):
                for otherBase2 in seqlib.allOtherBases(base2, rna):
                    variant = seq[:i] + otherBase1 + seq[i+1:i+j+1] + otherBase2 + seq[i+j+2:]
                    listVariants.append(variant)

    return listVariants

# Generate a list of variants with double bulges
def generate_1x0_1x0(consensus, rna=False):
    
    seq = seqlib.standardize(consensus, rna)
    listVariants = []
    for i, base1 in enumerate(seq[:-1]):
        for j, base2 in enumerate(seq[i+1:-1]):
            for newBase1 in seqlib.allBases(rna):
                for newBase2 in seqlib.allBases(rna):
                    variant = seq[:i+1] + newBase1 + seq[i+1:i+j+2] + newBase2 + seq[i+j+2:]
                    listVariants.append(variant)

    return listVariants

# Generate a list of variants with double deletions
def generate_0x1_0x1(consensus, rna=False):
    
    seq = seqlib.standardize(consensus, rna)
    listVariants = []
    for i, base1 in enumerate(seq):
        for j, base2 in enumerate(seq[i+1:]):
            variant = seq[:i] + seq[i+1:i+j+1] + seq[i+j+2:]
            listVariants.append(variant)

    return listVariants

###########################################################

# Check duplicity between two list of variants
def checkDuplicity(listVariants1, listVariants2):
    
    duplicates = 0
    
    for variant1 in listVariants1:
        for variant2 in listVariants2:
            if variant1 == variant2:
                duplicates += 1

    return duplicates
