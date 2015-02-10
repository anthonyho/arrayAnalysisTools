#!/usr/bin/env python
#
# Map barcodes
#
# Try to find the annotations of all barcodes in a given
# CPseq-type file from the barcode-annotation dictionary
#
# Anthony Ho, ahho@stanford.edu, 2/2/2015
# Last update 2/2/2015


# Import libraries
import os, sys
import argparse
import pandas as pd
import numpy as np


def main():

    # Get options and arguments from command line
    parser = argparse.ArgumentParser(description="map barcodes")
    parser.add_argument('bcDictFilePath', help="path to the barcode-annotation dictionary file")
    parser.add_argument('seqsFilePath', help="path to the CPseq-type file to be annotated")
    parser.add_argument('outputFilePath', help="path to the output file")
    args = parser.parse_args()

    # Read barcode-annotation file as Pandas dataframe
    bcDF = pd.read_csv(args.bcDictFilePath, sep='\t', usecols=['barcode', 'annotation'])

    # Check if any barcodes duplicate
    if np.sum(bcDF.duplicated('barcode')):
        print "There are barcode duplicates in the barcode-annotation dictionary!"
        print "Quitting now..."
        return 0

    # Convert barcode-annotation dataframe into dictionary
    bcDict = bcDF.set_index('barcode')['annotation'].to_dict()

    # Load seqFile
    seqsDF = pd.read_csv(args.seqsFilePath, sep='\t')

    # Map barcodes and represent missing keys entries with 'NA:nan:nan:nan:::'
    allAnnotations = seqsDF['barcode'].apply(bcDict.get, args=None).replace({None: 'NA:nan:nan:nan:::'})

    # Insert annotations to dataframe
    seqsDF.insert(0, 'annotation', allAnnotations)

    # Drop the clusterID column
    seqsDF = seqsDF.drop('clusterID', 1)

    # Write to output file path
    seqsDF.to_csv(args.outputFilePath, sep='\t', index=False)

    return 1

if __name__ == "__main__":
    main()
