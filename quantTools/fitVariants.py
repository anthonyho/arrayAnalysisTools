#!/usr/bin/env python
#
# Fit variants with curves
#
# Anthony Ho, ahho@stanford.edu, 11/29/2016
# Last update 9/2/2018
#

# Import libraries
import os, sys
import argparse
import numpy as np
import pandas as pd
import fitlib


def binding(params, conc):
    RT = 0.59248475334 # 25C
    dG = params[0]
    Kd = np.exp(dG / RT)
    fmin = params[1]
    fmax = params[2]
    return (fmax - fmin) / (1 + Kd / conc) + fmin

def jac_binding(params, conc):
    RT = 0.59248475334
    dG = params[0]
    Kd = np.exp(dG / RT)
    fmin = params[1]
    fmax = params[2]
    E = Kd / conc
    partial_dG = -E / (RT * (1 + E)**2)
    partial_b = - 1 / (1 + E) + 1
    partial_A = 1 / (1+E)
    return np.vstack([partial_dG, partial_b, partial_A]).T

def fit(row, conc, params0, bounds, method, sigma=None):
    return fitlib.lsqcurvefit(binding, conc, row, params0, sigma=sigma,
                              jac=jac_binding, bounds=bounds, method=method, disp=False)
    
def fit_all(conc, df_m,
            params0=[-1, 0, 0.8], bounds=[(None, -0.9), (0, 1), (0, 1)],
            method='TNC', df_s=None):
    list_fitObj = []
    for i in range(0, len(df_m)):
        if i % 10000 == 0:
            print "Fitting the {}th variant".format(i+1)
        if df_s is not None:
            fitObj = fit(df_m.iloc[i], conc, params0, bounds, method, sigma=df_s.iloc[i])
        else:
            fitObj = fit(df_m.iloc[i], conc, params0, bounds, method)
        list_fitObj.append(fitObj)
    return list_fitObj

def extract(list_fitObj, num):
    list_results = []
    for i, fitObj in enumerate(list_fitObj):
        curr_result = []
        curr_result += list(fitObj.params)
        curr_result += list(fitObj.paramSEs)
        curr_result += list(fitObj.paramPvals)
        curr_result.append(fitObj.R2)
        curr_result.append(fitObj.adjR2)
        curr_result.append(fitObj.reChi2)
        curr_result.append(fitObj.SER)
        curr_result.append(num['n'].iloc[i])
        list_results.append(curr_result)
    columns = ['dG', 'fmin', 'fmax',
               'dG_SE', 'fmin_SE', 'fmax_SE',
               'dG_pval', 'fmin_pval', 'fmax_SE_pval',
               'R2', 'adjR2', 'reChi2', 'SER', 'count']
    return pd.DataFrame(list_results, columns=columns, index=num.index)

def main():

    # Get options and arguments from command line
    parser = argparse.ArgumentParser(description="fit variants to curves")
    parser.add_argument('medCPvariantFilePath', help="path to the median CPvariant file")
    parser.add_argument('semCPvariantFilePath', help="path to the sem CPvariant file")
    parser.add_argument('countCPvariantFilePath', help="path to the count CPvariant file")
    parser.add_argument('concFilePath', help="path to the concentration file")
    parser.add_argument('outputFilePath', help="path to the concentration file")
    args = parser.parse_args()
    
    # Read CPvariant files
    median = pd.read_csv(args.medCPvariantFilePath,
                         sep='\t', index_col=0)
    sem = pd.read_csv(args.semCPvariantFilePath,
                      sep='\t', index_col=0)
    num = pd.read_csv(args.countCPvariantFilePath,
                      sep='\t', index_col=0)
    # Read conc
    conc = pd.read_csv(args.concFilePath).values.flatten()

    # Fit
    list_fitObj = fit_all(conc, median.iloc[:, 1:])
    
    # Extract
    df_results = extract(list_fitObj, num)

    # Save
    df_results.to_csv(args.outputFilePath, sep='\t')

    return 1

if __name__ == "__main__":
    main()
