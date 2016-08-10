#!/usr/bin/env python
#
# Fit single clusters
#
# Anthony Ho, ahho@stanford.edu, 2/3/2015
# Last update 2/24/2015


# Import libraries
import os, sys
import argparse
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
import multiprocessing
from itertools import izip
import fitlib
import parselib


def main():

    # Get options and arguments from command line
    parser = argparse.ArgumentParser(description="fit single clusters")
    parser.add_argument('-n', '--numCore', type=int, default=1, help="number of cores to use (default=1)")
    parser.add_argument('-v', '--verbose', type=int, default=0, help="verbosity of progress. 0 = no verbosity. (default=0)")
    parser.add_argument('fitParamFilePath', help="path to the file that specifies fitting parameters")
    parser.add_argument('inputFilePath', help="path to the file containing the raw signals of the clusters to be fitted")
    parser.add_argument('outputFilePath', help="path to the output file")
    args = parser.parse_args()

    # Define default attributes from the fit result to write to the output file
    outputAttrs = ['params',
                   'paramSEs',
                   'paramTvals',
                   'paramPvals',
                   'RSS',
                   'reChi2',
                   'SER',
                   'nit',
                   'status']

    # Import fit parameters from fitParamFile
    fitParamDict = fitlib.lsqcurvefit.parseFitParamFromFile(args.fitParamFilePath)

    # Read inputFile as Pandas dataframe
    allClusters = pd.read_csv(args.inputFilePath, sep='\t')
    allSignals = parselib.splitConcatedDFColumnIntoNDarray(allClusters['signals'], ':')
    allTimes = parselib.splitConcatedDFColumnIntoNDarray(allClusters['times'], ':')

    # Fit single clusters
    fitResults = Parallel(n_jobs=args.numCore, verbose=args.verbose)(delayed(fitlib.lsqcurvefit)(x=time, y=signal, **fitParamDict)
                                                                     for time, signal in izip(allTimes, allSignals))

    # Add attributes as defined in outputAttrs as columns in the allClusters dataframe
    for attr in outputAttrs:
        val0 = getattr(fitResults[0], attr)
        if not isinstance(val0, np.ndarray):
            allClusters[attr] = [getattr(r, attr) for r in fitResults]
        else:
            attrNames = [attr+str(i+1) for i in range(len(val0))]
            allClusters = allClusters.join(pd.DataFrame([getattr(r, attr) for r in fitResults], columns=attrNames))

    allClusters.to_csv(args.outputFilePath, sep='\t', index=False)

    return 1

if __name__ == "__main__":
    main()
