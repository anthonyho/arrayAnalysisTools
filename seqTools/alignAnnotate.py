#!/usr/bin/env python
#
# Align and annotate a given list of sequences
# against a set of reference sequences
#
# Anthony Ho, ahho@stanford.edu, 1/28/2015
# Last update 1/28/2015


## Import libraries
import os, sys
import argparse
import pandas as pd
import numpy as np
import alignLibs


def main():

    ## Get options and arguments from command line
    parser = argparse.ArgumentParser(description="align and annotate a given list of sequences against a set of reference sequences")
    parser.add_argument('-i', '--indel', action='store_true', help="enable alignment and annotation of indel sequences")
    parser.add_argument('refSeqFilePath', help="path to the reference sequences file (in FASTQ format)")
    parser.add_argument('seqCol', help="column of the file containing the list of sequences to be aligned and annotated (in Python notation)")
    parser.add_argument('seqFilePath', help="path to the file containing the list of sequences to be aligned and annotated")
    parser.add_argument('outCol', help="column of the output file to which the annotation will be written to")
    parser.add_argument('outputFilePath', help="path to the output file")
    args = parser.parse_args()


    ## Initialization

    


    




    
    return 1 

if __name__ == "__main__":
    main()
