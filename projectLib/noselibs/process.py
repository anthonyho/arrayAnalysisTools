# Anthony Ho, ahho@stanford.edu, 1/5/2016
# Last update 2/6/2015
"""Library containing functions to process data"""


import numpy as np
import pandas as pd
import aux


# Merge and process all variants related information into a single multilabel dataframe
def mergeAllVariants(variants_dict, variants_max_dict, bindingSeries_dict, concentrations_dict, annotatedClusters):

    # Merge variants_dict and variants_max_dict into multilabel dataframes
    variants = pd.concat(variants_dict, axis=1)
    variants_max = pd.concat({ligand: variants_max_dict[ligand].rename(columns={'0': 'max'}) for ligand in variants_max_dict}, axis=1)
    
    # Compute median signals of binding series for each variant, and merge dict of binding series into a single multilabel dataframe
    medianBindingSeriesByVariants_dict = {}
    medianNormBindingSeriesByVariants_dict = {}
    ciBindingSeriesByVariants_dict = {}
    ciNormBindingSeriesByVariants_dict = {}
    for ligand in bindingSeries_dict:
        
        groupedBindingSeries = pd.merge(bindingSeries_dict[ligand], annotatedClusters,
                                        how='inner', left_index=True, right_index=True).groupby('variant_number')
        
        medianBindingSeriesByVariants = groupedBindingSeries.median()
        medianBindingSeriesByVariants.columns = ['bs_{}'.format(c) for c in concentrations_dict[ligand]/1000]
        medianBindingSeriesByVariants_dict[ligand] = medianBindingSeriesByVariants

        medianNormBindingSeriesByVariants = medianBindingSeriesByVariants.divide(variants_max[ligand]['max'], axis=0)
        medianNormBindingSeriesByVariants.columns = [col+'_norm' for col in medianNormBindingSeriesByVariants.columns]
        medianNormBindingSeriesByVariants_dict[ligand] = medianNormBindingSeriesByVariants

        ciBindingSeriesByVariants = groupedBindingSeries.std() / np.sqrt(groupedBindingSeries.count())
        ciBindingSeriesByVariants.columns = ['bs_{}_err'.format(c) for c in concentrations_dict[ligand]/1000]
        ciBindingSeriesByVariants_dict[ligand] = ciBindingSeriesByVariants

        ciNormBindingSeriesByVariants = ciBindingSeriesByVariants.divide(variants_max[ligand]['max'], axis=0)
        ciNormBindingSeriesByVariants.columns = [col+'_norm' for col in ciNormBindingSeriesByVariants.columns]
        ciNormBindingSeriesByVariants_dict[ligand] = ciNormBindingSeriesByVariants
    
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

    computedColumns = {'Kd': aux.dGtoKd(variants_all_swapped['dG'], unit='uM'),
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

    mixtures_all = pd.concat({cm: mixtures_dict[cm][[0, 2]].rename(columns={'0': 'cs_norm', '2': 'cs'}) for cm in mixtures_dict}, axis=1)

    return mixtures_all


# Filter variants
def filterVariants(df, dG_max=None, dG_err_max=None, 
                   fmax_max=None, fmax_norm_max=None, fmax_err_max=None, fmax_err_norm_max=None, frac_fmax_err_max=None,
                   rsq_min=None, nCluster_min=None, pvalue_max=None, minColsAgreed=1, swapped=False):

    if minColsAgreed == 'all':
        minColsAgreed = len(df.columns.get_level_values(0).unique())
    
    if isinstance(df.columns, pd.core.index.MultiIndex):
        df_swapped = df.swaplevel(0, 1, axis=1).sort_index(axis=1)
    else:
        df_swapped = df

    filters = df_swapped['dG'] > -np.Inf
    if dG_max is not None:
        filters = filters & (df_swapped['dG'] < dG_max)
    if dG_err_max is not None:
        filters = filters & (df_swapped['dG_err'] < dG_err_max)
    if fmax_max is not None:
        filters = filters & (df_swapped['fmax'] < fmax_max)
    if fmax_norm_max is not None:
        filters = filters & (df_swapped['fmax_norm'] < fmax_norm_max)
    if fmax_err_max is not None:
        filters = filters & (df_swapped['fmax_err'] < fmax_err_max)
    if fmax_err_norm_max is not None:
        filters = filters & (df_swapped['fmax_err_norm'] < fmax_err_norm_max)
    if frac_fmax_err_max is not None:
        filters = filters & ((df_swapped['fmax_err'] / df_swapped['fmax']) < frac_fmax_err_max)
    if rsq_min is not None:
        filters = filters & (df_swapped['rsq'] > rsq_min)
    if nCluster_min is not None:
        filters = filters & (df_swapped['numClusters'] >= nCluster_min)
    if pvalue_max is not None:
        filters = filters & (df_swapped['pvalue'] < pvalue_max)
        
    filtered_index = filters.sum(axis=1) >= minColsAgreed

    if swapped:
        return df_swapped[filtered_index]
    else:
        return df_swapped[filtered_index].swaplevel(0, 1, axis=1).sort_index(axis=1)
