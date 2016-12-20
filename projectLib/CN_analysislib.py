# Anthony Ho, ahho@stanford.edu, 8/31/2016
# Last update 10/11/2016
"""Python module containing analysis functions for chemical nose project"""


import numpy as np
import pandas as pd
import varlib
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis


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
def filterVariants(df, dG_max=None, dG_err_max=None, rsq_min=None, nCluster_min=None, pvalue_max=None, minColsAgreed=1):

    if minColsAgreed == 'all':
        minColsAgreed = len(df.columns.get_level_values(0).unique())
    
    if isinstance(df.columns, pd.core.index.MultiIndex):
        df_swapped = df.swaplevel(0, 1, axis=1).sort_index(axis=1)
    else:
        df_swapped = df

    filters = pd.Series([True] * len(df_swapped), index=df_swapped.index)
    if dG_max is not None:
        filters = filters & ((df_swapped['dG'] < dG_max).sum(axis=1) >= minColsAgreed)
    if dG_err_max is not None:
        filters = filters & (((df_swapped['dG_ub'] - df_swapped['dG_lb']) < dG_err_max).sum(axis=1) >= minColsAgreed)
    if rsq_min is not None:
        filters = filters & ((df_swapped['rsq'] > rsq_min).sum(axis=1) >= minColsAgreed)
    if nCluster_min is not None:
        filters = filters & ((df_swapped['numClusters'] >= nCluster_min).sum(axis=1) >= minColsAgreed)
    if pvalue_max is not None:
        filters = filters & ((df_swapped['pvalue'] < pvalue_max).sum(axis=1) >= minColsAgreed)

    return df_swapped[filters].swaplevel(0, 1, axis=1).sort_index(axis=1)


# Perform PCA with sklearn's PCA to reduce dimensionality of aptamers
def performPCA(data_to_fit, numComponent=None):
    data_to_fit_np_t = np.array(data_to_fit).T
    if numComponent is None:
        numComponent = len(data_to_fit_np_t)    
    pca_model = PCA(n_components=numComponent, whiten=False)
    pca_results = pca_model.fit_transform(data_to_fit_np_t)
    return pca_model, pca_results


# Perform LDA with sklearn's LDA to reduce dimensionality of aptamers
def performLDA(data_to_fit, y, numComponent=None):
    data_to_fit_np_t = np.array(data_to_fit).T
    if numComponent is None:
        numComponent = len(data_to_fit_np_t)
    lda_model = LinearDiscriminantAnalysis(n_components=numComponent)
    lda_results = lda_model.fit_transform(data_to_fit_np_t, y)
    return lda_model, lda_results


# Perform LDA with sklearn's LDA to reduce dimensionality of aptamers
def performPCA_LDA(data_to_fit, y, numPCAComponent=None, numLDAComponent=None):

    pca_model, pca_results = performPCA(data_to_fit, numComponent=numPCAComponent)
    lda_model, lda_results = performLDA(pca_results.T, y, numComponent=numLDAComponent)

    return pca_model, pca_results, lda_model, lda_results


