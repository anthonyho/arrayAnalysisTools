# Anthony Ho, ahho@stanford.edu, 1/5/2016
# Last update 1/5/2015
"""Library containing functions to process data"""


import numpy as np
import pandas as pd
import liblib


# Merge and process all variants related information into a single multilabel dataframe
def mergeAllVariants(variants_dict, variants_max_dict, bindingSeries_dict, concentrations_dict, annotatedClusters):

    # Merge variants_dict and variants_max_dict into multilabel dataframes
    variants = pd.concat(variants_dict, axis=1)
    variants_max = pd.concat({currSM: variants_max_dict[currSM].rename(columns={'0': 'max'}) for currSM in variants_max_dict}, axis=1)
    
    # Compute median signals of binding series for each variant, and merge dict of binding series into a single multilabel dataframe
    medianBindingSeriesByVariants_dict = {}
    medianNormBindingSeriesByVariants_dict = {}
    ciBindingSeriesByVariants_dict = {}
    ciNormBindingSeriesByVariants_dict = {}
    for currSM in bindingSeries_dict:
        
        groupedBindingSeries = pd.merge(bindingSeries_dict[currSM], annotatedClusters,
                                        how='inner', left_index=True, right_index=True).groupby('variant_number')
        
        medianBindingSeriesByVariants = groupedBindingSeries.median()
        medianBindingSeriesByVariants.columns = ['bs_'+str(c) for c in concentrations_dict[currSM]/1000]
        medianBindingSeriesByVariants_dict[currSM] = medianBindingSeriesByVariants

        medianNormBindingSeriesByVariants = medianBindingSeriesByVariants.divide(variants_max[currSM]['max'], axis=0)
        medianNormBindingSeriesByVariants.columns = [col+'_norm' for col in medianNormBindingSeriesByVariants.columns]
        medianNormBindingSeriesByVariants_dict[currSM] = medianNormBindingSeriesByVariants

        ciBindingSeriesByVariants = groupedBindingSeries.std() / np.sqrt(groupedBindingSeries.count())
        ciBindingSeriesByVariants.columns = ['ci_bs_'+str(c) for c in concentrations_dict[currSM]/1000]
        ciBindingSeriesByVariants_dict[currSM] = ciBindingSeriesByVariants

        ciNormBindingSeriesByVariants = ciBindingSeriesByVariants.divide(variants_max[currSM]['max'], axis=0)
        ciNormBindingSeriesByVariants.columns = [col+'_norm' for col in ciNormBindingSeriesByVariants.columns]
        ciNormBindingSeriesByVariants_dict[currSM] = ciNormBindingSeriesByVariants
    
    variants_bindingSeries = pd.concat(medianBindingSeriesByVariants_dict, axis=1)
    variants_normBindingSeries = pd.concat(medianNormBindingSeriesByVariants_dict, axis=1)
    variants_ciBindingSeries = pd.concat(ciBindingSeriesByVariants_dict, axis=1)
    variants_ciNormBindingSeries = pd.concat(ciNormBindingSeriesByVariants_dict, axis=1)

    # Merge variants_dict, variants_max_dict, and bindingSeries into a single multilabel dataframe
    variants_all = pd.concat([variants, variants_max, 
                              variants_bindingSeries, variants_normBindingSeries, 
                              variants_ciBindingSeries, variants_ciNormBindingSeries], axis=1).sort_index(axis=1)

    # Compute new columns
    variants_all_swapped = variants_all.swaplevel(0, 1, axis=1).sort_index(axis=1)

    computedColumns = {'Kd': liblib.dGtoKd(variants_all_swapped['dG'], unit='uM'),
                       'dG_err': (variants_all_swapped['dG_ub'] - variants_all_swapped['dG_lb']) / 3.92,
                       'fmax_err': (variants_all_swapped['fmax_ub'] - variants_all_swapped['fmax_lb']) / 3.92,
                       'fmax_norm': variants_all_swapped['fmax'] / variants_all_swapped['max'],
                       'fmax_err_norm': (variants_all_swapped['fmax_ub'] - variants_all_swapped['fmax_lb']) / 3.92 / variants_all_swapped['max']
                       }
    
    df_computedColumns = pd.concat(computedColumns, axis=1).swaplevel(0, 1, axis=1).sort_index(axis=1)
    
    variants_all = pd.concat([variants_all, df_computedColumns], axis=1).sort_index(axis=1)
    
    return variants_all


# Merge mixtures_dict into a single multilabel dataframe
def mergeAllMixtures(mixtures_dict):

    mixtures_all = pd.concat({currCM: mixtures_dict[currCM][[0, 2]].rename(columns={'0': 'N', '2': 'S'}) for currCM in mixtures_dict}, axis=1)

    return mixtures_all


# Filter variants
def filterVariants(df, dG_max=None, dG_err_max=None, rsq_min=None, nCluster_min=None, pvalue_max=None, minColsAgreed=1, swapped=False):

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

    if swapped:
        return df_swapped[filters]
    else:
        return df_swapped[filters].swaplevel(0, 1, axis=1).sort_index(axis=1)

