# Anthony Ho, ahho@stanford.edu, 8/31/2016
# Last update 8/31/2016
"""Python module containing analysis functions for chemical nose project"""


import numpy as np
import pandas as pd
import varlib


# Generate library sequences
def generateLibSeq(sublib='all', form='dict'):
    
    # Define names and sequences
    LP5a_refSeq = "NNGGATTTTCCNNNACGAAGTGTCCCGAGA"
    LP5b_refSeq = "NNGGATTTTCCNNACGAAGTNGTCCCGAGA"
    LP6_refSeq =  "NNGGATTTTCCNNNACGAAGTNGTCCCGAG"
    LP7a_refSeq = "NNGGATTTTCCNNNNACGAAGTNGTCCCGA"
    LP7b_refSeq = "NNGGATTTTCCNNNACGAAGTNNGTCCCGA"
    LP8_refSeq =  "NNGGATTTTCCNNNNACGAAGTNNGTCCCG"
    allRefSeqs = [LP5a_refSeq, LP5b_refSeq, LP6_refSeq, LP7a_refSeq, LP7b_refSeq, LP8_refSeq]
    libNames = ['LP5a', 'LP5b', 'LP6', 'LP7a', 'LP7b', 'LP8']

    # Generate library sequences
    libSeqs_dict = {}
    for name, refSeq in zip(libNames, allRefSeqs):
        libSeqs_dict[name] = varlib.generateSatLib(refSeq)

    # Return the appropiate sub library or form
    if sublib != 'all':
        return libSeqs_dict[sublib]
    else:
        if form == 'dict':
            return libSeqs_dict
        elif form == 'list':
            return [seq for name in libNames for seq in libSeqs_dict[name]]
        elif form == 'df':
            listAllLibSeqs = [seq for name in libNames for seq in libSeqs_dict[name]]
            listAllLibNames = []
            for name in libNames:
                listAllLibNames = listAllLibNames + [name] * len(libSeqs_dict[name])
            return pd.DataFrame({'lib': listAllLibNames, 'seq': listAllLibSeqs})


# Filter variants
def filterVariants(data, dG_err_max=None, rsq_min=None, nCluster_min=None, pvalue_max=None):
    
    filters = pd.Series([True] * len(data), index=data.index)
    if dG_err_max is not None:
        filters = filters & ((data['dG_ub'] - data['dG_lb']) < dG_err_max)
    if rsq_min is not None:
        filters = filters & (data['rsq'] > rsq_min)
    if nCluster_min is not None:
        filters = filters & (data['numClusters'] >= nCluster_min)
    if pvalue_max is not None:
        filters = filters & (data['pvalue'] < pvalue_max)

    return data[filters]
