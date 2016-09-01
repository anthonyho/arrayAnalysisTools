# Anthony Ho, ahho@stanford.edu, 8/31/2016
# Last update 8/31/2016
"""Python module containing analysis functions for chemical nose project"""


import numpy as np
import pandas as pd


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
