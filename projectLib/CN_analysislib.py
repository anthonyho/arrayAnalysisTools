# Anthony Ho, ahho@stanford.edu, 8/31/2016
# Last update 8/31/2016
"""Python module containing analysis functions for chemical nose project"""


import numpy as np
import pandas as pd
import varlib


# Generate library sequences
def generateLibSeq(sublib='all', form='dict'):
    
    # Define names and reference sequences
    libNames = ['LP5a', 'LP5b', 'LP6', 'LP7a', 'LP7b', 'LP8']
    refSeqs = {'LP5a': "NNGGATTTTCCNNNACGAAGTGTCCCGAGA", 
               'LP5b': "NNGGATTTTCCNNACGAAGTNGTCCCGAGA", 
               'LP6':  "NNGGATTTTCCNNNACGAAGTNGTCCCGAG", 
               'LP7a': "NNGGATTTTCCNNNNACGAAGTNGTCCCGA", 
               'LP7b': "NNGGATTTTCCNNNACGAAGTNNGTCCCGA", 
               'LP8':  "NNGGATTTTCCNNNNACGAAGTNNGTCCCG"}

    # Generate library sequences
    libSeqs_dict = {}
    for name in libNames:
        libSeqs_dict[name] = varlib.generateSatLib(refSeqs[name])

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


# Return the sublibrary of which a library sequence belongs to
def returnSubLib(seq):

    # Get full list of library sequences with their sub-library
    libSeqs = generateLibSeq(sublib='all', form='df')
    name = libSeqs.set_index('seq').loc[seq]['lib']
    
    return name


# Return the reference sequence given a library sequence
def returnRefSeq(seq):

    # Define reference sequences
    refSeqs = {'LP5a': "NNGGATTTTCCNNNACGAAGTGTCCCGAGA", 
               'LP5b': "NNGGATTTTCCNNACGAAGTNGTCCCGAGA", 
               'LP6':  "NNGGATTTTCCNNNACGAAGTNGTCCCGAG", 
               'LP7a': "NNGGATTTTCCNNNNACGAAGTNGTCCCGA", 
               'LP7b': "NNGGATTTTCCNNNACGAAGTNNGTCCCGA", 
               'LP8':  "NNGGATTTTCCNNNNACGAAGTNNGTCCCG"}

    # Get full list of library sequences with their sub-library
    libSeqs = generateLibSeq(sublib='all', form='df')
    name = libSeqs.set_index('seq').loc[seq]['lib']
    
    return refSeqs[name]


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
