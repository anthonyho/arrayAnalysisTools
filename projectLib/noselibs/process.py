# Anthony Ho, ahho@stanford.edu, 1/5/2017
# Last update 3/15/2017
"""Library containing functions to process data"""


from __future__ import division
import numpy as np
import pandas as pd
import aux


# Conversion factor of one sided 95% CI to standard error
cf = 1.96


# Merge and process all variants related information
# into a single multilabel dataframe
def merge_variants(variants_dict, variants_max_dict, binding_series_dict,
                   conc_dict, annotated_clusters):

    # Merge variants_dict and variants_max_dict into multilabel dataframes
    variants = pd.concat(variants_dict, axis=1)
    variants_max = pd.concat({ligand: variants_max_dict[ligand].rename(columns={'0': 'max'})
                              for ligand in variants_max_dict}, axis=1)

    # Compute median signals of binding series for each variant,
    # and merge dict of binding series into a single multilabel dataframe
    med_bs_by_var_dict = {}
    med_norm_bs_by_var_dict = {}
    err_bs_by_var_dict = {}
    err_norm_bs_by_var_dict = {}
    for ligand in binding_series_dict:

        bs = binding_series_dict[ligand]
        conc = conc_dict[ligand] / 1000
        vmax = variants_max[ligand]['max']

        bs_by_var = pd.merge(bs, annotated_clusters,
                             left_index=True, right_index=True,
                             how='inner').groupby('variant_number')

        med_bs_by_var = bs_by_var.median()
        med_bs_by_var.columns = ['bs_{}'.format(c) for c in conc]
        med_bs_by_var_dict[ligand] = med_bs_by_var

        med_norm_bs_by_var = med_bs_by_var.divide(vmax, axis=0)
        med_norm_bs_by_var.columns = [col+'_norm'
                                      for col in med_norm_bs_by_var.columns]
        med_norm_bs_by_var_dict[ligand] = med_norm_bs_by_var

        err_bs_by_var = bs_by_var.std() / np.sqrt(bs_by_var.count())
        err_bs_by_var.columns = ['bs_{}_err'.format(c) for c in conc]
        err_bs_by_var_dict[ligand] = err_bs_by_var

        err_norm_bs_by_var = err_bs_by_var.divide(vmax, axis=0)
        err_norm_bs_by_var.columns = [col+'_norm'
                                      for col in err_norm_bs_by_var.columns]
        err_norm_bs_by_var_dict[ligand] = err_norm_bs_by_var

    variants_bs = pd.concat(med_bs_by_var_dict, axis=1)
    variants_norm_bs = pd.concat(med_norm_bs_by_var_dict, axis=1)
    variants_err_bs = pd.concat(err_bs_by_var_dict, axis=1)
    variants_err_norm_bs = pd.concat(err_norm_bs_by_var_dict, axis=1)

    # Merge variants_dict, variants_max_dict, and bindingSeries
    # into a single multilabel dataframe
    variants_all = pd.concat([variants, variants_max,
                              variants_bs, variants_norm_bs,
                              variants_err_bs, variants_err_norm_bs],
                             axis=1).sort_index(axis=1)

    # Compute new columns
    swapped = variants_all.swaplevel(0, 1, axis=1).sort_index(axis=1)

    imputed = {'Kd': aux.dG_to_Kd(swapped['dG'], unit='uM'),
               'dG_err': (swapped['dG'] - swapped['dG_lb']) / cf,
               'fmax_err': (swapped['fmax'] - swapped['fmax_lb']) / cf,
               'fmax_norm': swapped['fmax'] / swapped['max'],
               'fmax_err_norm': (swapped['fmax'] - swapped['fmax_lb']) / swapped['max'] / cf,
               'fmax_err_frac': (swapped['fmax'] - swapped['fmax_lb']) / swapped['fmax'] / cf
               }

    imputed = pd.concat(imputed, axis=1)

    variants_all = pd.concat([swapped, imputed],
                             axis=1).swaplevel(0, 1, axis=1).sort_index(axis=1)

    return variants_all


# Merge mixtures_dict into a single multilabel dataframe
def merge_mixtures(mixtures_dict,
                   cols=[0, 2], names={'0': 'cs_norm', '2': 'cs'}):

    mixtures_all = pd.concat({cm: mixtures_dict[cm][cols].rename(columns=names)
                              for cm in mixtures_dict}, axis=1)

    return mixtures_all


# Filter variants
def filter_variants(df, dG_max=None, dG_err_max=None,
                    fmax_max=None, fmax_norm_max=None,
                    fmax_err_max=None, fmax_err_norm_max=None,
                    frac_fmax_err_max=None,
                    rsq_min=None, n_cluster_min=None, pvalue_max=None,
                    min_cols_agreed=1, swapped=False):

    if min_cols_agreed == 'all':
        min_cols_agreed = len(df.columns.get_level_values(0).unique())

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
        filters = filters & (df_swapped['fmax_err_frac'] < frac_fmax_err_max)
    if rsq_min is not None:
        filters = filters & (df_swapped['rsq'] > rsq_min)
    if n_cluster_min is not None:
        filters = filters & (df_swapped['numClusters'] >= n_cluster_min)
    if pvalue_max is not None:
        filters = filters & (df_swapped['pvalue'] < pvalue_max)

    filtered_index = filters.sum(axis=1) >= min_cols_agreed

    if swapped:
        return df_swapped[filtered_index]
    else:
        return df_swapped[filtered_index].swaplevel(0, 1, axis=1).sort_index(axis=1)
