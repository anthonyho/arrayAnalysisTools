#!/usr/bin/env python
#
# Combine CPseq-type tiles into a single file and add
# a column indicating which tile the cluster is from
#
# Anthony Ho, ahho@stanford.edu, 2/26/2015
# Last update 2/26/2015


# Import libraries
import os, sys
import argparse
import pandas as pd
import numpy as np


def main():

    # Get options and arguments from command line
    parser = argparse.ArgumentParser(description="combine tiles")
    parser.add_argument('outputFilePath', help="path to the output file")
    parser.add_argument('tileFilePaths', nargs='*', help="paths to the tile files to be combined")
    args = parser.parse_args()

    # Import every tile into separate dataframe
    listTilesDF = []

    for tile in args.tileFilePaths:

        # Get tile number
        (_, tileFilename) = os.path.split(tile)
        (tileBasename, _) = os.path.splitext(tileFilename)
        tileNumber = int(tileBasename.split('tile')[1][0:3])

        # Load dataframe
        tileDF = pd.read_csv(tile, sep='\t')

        # Insert tile column into dataframe
        tileDF.insert(2, 'tile', [tileNumber]*len(tileDF.index))

        listTilesDF.append(tileDF)

    # Concat into a single dataframe
    allTilesDF = pd.concat(listTilesDF)

    # Write to output file path
    allTilesDF.to_csv(args.outputFilePath, sep='\t', index=False)

    return 1

if __name__ == "__main__":
    main()
