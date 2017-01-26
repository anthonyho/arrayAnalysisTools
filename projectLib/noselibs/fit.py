# Anthony Ho, ahho@stanford.edu, 1/5/2017
# Last update 1/25/2017
"""Library containing the switching equations, their derivatives, residuals, and fitting functions"""


import numpy as np
import pandas as pd
import lmfit
import liblib


# --- Private library functions --- #

RT = 0.582


def _switchingEq(mu, dG, fmax, fmin):
    '''Compute the intensity of a switching aptamer'''
    B = np.exp(-(dG - mu)/RT)
    Q = B.sum(axis=1).reshape(-1, 1)
    Q_weighted = (B * (fmax - fmin)).sum(axis=1).reshape(-1, 1)
    return (Q_weighted / (1 + Q) + fmin).reshape(-1)


def _switchingEq_jaconbian_mu(mu, dG, fmax, fmin):
    '''Compute the Jacobian of the switching equation with respect to conc'''
    B = np.exp(-(dG - mu)/RT)
    Q = B.sum(axis=1).reshape(-1, 1)
    Q_weighted = (B * (fmax - fmin)).sum(axis=1).reshape(-1, 1)
    return B / RT * ((fmax - fmin) * (1 + Q) - Q_weighted) / (1 + Q)**2


def _switchingEq_jacobian_dG(mu, dG, fmax, fmin):
    '''Compute the Jacobian of the switching equation with respect to dG'''
    B = np.exp(-(dG - mu)/RT)
    Q = B.sum(axis=1).reshape(-1, 1)
    Q_weighted = (B * (fmax - fmin)).sum(axis=1).reshape(-1, 1)
    return - B / RT * ((fmax - fmin) * (1 + Q) - Q_weighted) / (1 + Q)**2


def _switchingEq_jacobian_fmax(mu, dG):
    '''Compute the Jacobian of the switching equation with respect to fmax'''
    B = np.exp(-(dG - mu)/RT)
    Q = B.sum(axis=1).reshape(-1, 1)
    return B / (1 + Q)


def _switchingEq_jacobian_fmin(mu, dG):
    '''Compute the Jacobian of the switching equation with respect to fmin'''
    B = np.exp(-(dG - mu)/RT)
    Q = B.sum(axis=1).reshape(-1, 1)
    return 1 / (1 + Q)


def _switchingEq_errors(mu, dG, fmax, fmin, data_err, dG_err, fmax_err, fmin_err):
    '''Compute the total uncertainty of each aptamer'''
    var_data = (data_err**2).reshape(-1)
    var_dG_sum = ((_switchingEq_jacobian_dG(mu, dG, fmax, fmin) * dG_err)**2).sum(axis=1)
    var_fmax_sum = ((_switchingEq_jacobian_fmax(mu, dG) * fmax_err)**2).sum(axis=1)
    var_fmin = ((_switchingEq_jacobian_fmin(mu, dG) * fmin_err)**2).reshape(-1)
    return np.sqrt(var_data + var_dG_sum + var_fmax_sum + var_fmin)


def _switchingEq_residuals(params, data, dG, fmax, fmin,
                           data_err, dG_err=None, fmax_err=None, fmin_err=None):
    '''Compute the residuals of each aptamer'''
    # Extract concentrations of ligands from params dict
    A, mu = _extractParams(params)

    # Reture residuals of each aptamers 
    if dG_err is None:
        return data - A * _switchingEq(mu, dG, fmax, fmin)
    else:
        return (data - A * _switchingEq(mu, dG, fmax, fmin)) / _switchingEq_errors(mu, dG, fmax, fmin, data_err, dG_err, fmax_err, fmin_err)


def _extractParams(params):
    '''Extract concentrations of ligands from params dict'''
    if isinstance(params, lmfit.Parameters):
        paramvals = params.valuesdict()
        A = paramvals.pop('A')
        mu = np.zeros((1, len(paramvals)))
        for i, currSM in enumerate(paramvals):
            mu[0, i] = paramvals[currSM]
    else:
        A = params[0]
        mu = np.array(params[1:]).reshape(1, -1)
    return A, mu


def _prepareVariables(data, dG, fmax=None, fmin=None,
                      data_err=None, dG_err=None, fmax_err=None, fmin_err=None):
    '''Initialize and typecast variables'''
    # Define fmax, fmin, their respective errors, and data_err if not provided
    if fmax is None:
        fmax = np.ones(dG.shape)
    if fmin is None:
        fmin = np.zeros(len(dG))
    if data_err is None:
        data_err = np.zeros(data.shape)
    if fmax_err is None:
        fmax_err = np.zeros(dG.shape)
    if fmin_err is None:
        fmin_err = np.zeros(len(dG))

    # Typecast and reshape input variables
    _data = np.array(data)
    _dG = np.array(dG)
    _fmax = np.array(fmax)
    _fmin = np.array(fmin).reshape(-1, 1)
    _data_err = np.array(data_err)
    _dG_err = np.array(dG_err)
    _fmax_err = np.array(fmax_err)
    _fmin_err = np.array(fmin_err).reshape(-1, 1)

    return _data, _dG, _fmax, _fmin, _data_err, _dG_err, _fmax_err, _fmin_err


# --- Functions meant for one sample --- #


def deconvoluteMixtures(data, dG, fmax=None, fmin=None, 
                        data_err=None, dG_err=None, fmax_err=None, fmin_err=None,
                        varyA=False, conc_init=None, conc_init_quantile=0.1, unit='uM', 
                        maxfev=500000, **kwargs):
    '''Fit the concentrations of ligands using lmfit'''
    # Initialize and typecast variables
    _data, _dG, _fmax, _fmin, _data_err, _dG_err, _fmax_err, _fmin_err = _prepareVariables(data, dG, fmax, fmin,
                                                                                           data_err, dG_err, fmax_err, fmin_err)

    # Define concenrations
    if conc_init is None or conc_init == 'auto':
        conc_init = np.array(liblib.dGtoKd(dG, unit=unit).quantile(conc_init_quantile))
    elif isinstance(conc_init, int) or isinstance(conc_init, float):
        conc_init = np.ones(len(dG.columns)) * conc_init
    mu_init = liblib.KdtodG(conc_init, unit=unit)
    params = lmfit.Parameters()
    params.add('A', value=1.0, min=0.0, vary=varyA)
    for i, currSM in enumerate(dG.columns):
        params.add(currSM, value=mu_init[i])
    
    # Fit and extract params
    fitResult = lmfit.minimize(_switchingEq_residuals, params,
                               args=(_data, _dG, _fmax, _fmin, _data_err, _dG_err, _fmax_err, _fmin_err), 
                               maxfev=maxfev, **kwargs)
    predictedConc = liblib.dGtoKd(pd.Series(fitResult.params.valuesdict()).drop('A'), unit=unit)
    
    return fitResult, predictedConc


def reportFit(output, weighted, 
              params, data, dG, fmax=None, fmin=None, 
              data_err=None, dG_err=None, fmax_err=None, fmin_err=None):
    '''Return various functional terms using the fitted parameters'''
    # Initialize and typecast variables
    _data, _dG, _fmax, _fmin, _data_err, _dG_err, _fmax_err, _fmin_err = _prepareVariables(data, dG, fmax, fmin,
                                                                                           data_err, dG_err, fmax_err, fmin_err)

    # Extract concentrations of ligands from params dict
    A, mu = _extractParams(params)
    
    # Output
    if output=='data':
        if weighted:
            return _data / _switchingEq_errors(mu, _dG, _fmax, _fmin, _data_err, _dG_err, _fmax_err, _fmin_err)
        else:
            return _data
    elif output=='eq':
        if weighted:
            return A * _switchingEq(mu, _dG, _fmax, _fmin) / _switchingEq_errors(mu, _dG, _fmax, _fmin, _data_err, _dG_err, _fmax_err, _fmin_err)
        else:
            return A * _switchingEq(mu, _dG, _fmax, _fmin)
    elif output=='residual':
        if weighted:
            return (_data - A * _switchingEq(mu, _dG, _fmax, _fmin)) / _switchingEq_errors(mu, _dG, _fmax, _fmin, _data_err, _dG_err, _fmax_err, _fmin_err)
        else:
            return _data - A * _switchingEq(mu, _dG, _fmax, _fmin)
    elif output=='err':
        return _switchingEq_errors(mu, _dG, _fmax, _fmin, _data_err, _dG_err, _fmax_err, _fmin_err)
    elif output=='err_mu':
        return _switchingEq_jacobian_mu(mu, _dG, _fmax, _fmin)
    elif output=='err_dG':
        return _switchingEq_jacobian_dG(mu, _dG, _fmax, _fmin)
    elif output=='err_fmax':
        return _switchingEq_jacobian_fmax(mu, _dG)
    elif output=='err_fmin':
        return _switchingEq_jacobian_fmin(mu, _dG)


# --- Higher level wrapper functions for multiple samples --- # 


def fitAllPureSamples(variants_subset, currConc, listSM, fmax=True, fmin=True, data_err=True, norm=True, 
                      varyA=False, conc_init=None, conc_init_quantile=0.1, **kwargs):
    '''Fit all the pure sample measurements at a given concentration'''
    # Define column names for fmax, fmin, and bs in the normalized and unnormalzied cases
    if norm:
        fmax_key = 'fmax_norm'
        fmax_err_key = 'fmax_err_norm'
        fmin_key = 'fmin_norm'
        fmin_err_key = 'fmin_err_norm'
        bs_key = 'bs_'+str(currConc)+'.0_norm'
        ci_bs_key = 'ci_bs_'+str(currConc)+'.0_norm'
    else:
        fmax_key = 'fmax'
        fmax_err_key = 'fmax_err'
        fmin_key = 'fmin'
        fmin_err_key = 'fmin_err'
        bs_key = 'bs_'+str(currConc)+'.0'
        ci_bs_key = 'ci_bs_'+str(currConc)+'.0'
    
    # Get fmax and fmin if requested
    if fmax:
        fmax = variants_subset[fmax_key]
        fmax_err = variants_subset[fmax_err_key]
    else:
        fmax = None
        fmax_err = None
    
    if fmin:
        fmin = variants_subset[fmin_key]
        fmin_err = variants_subset[fmin_err_key]
    else:
        fmin = None
        fmin_err = None
    
    # Fit all pure samples
    list_predictedConc = []
    dict_fitResults = {}
    for currSM in listSM:
        if data_err:
            fitResult, predictedConc = deconvoluteMixtures(variants_subset[bs_key][currSM], variants_subset['dG'], fmax, fmin,
                                                           variants_subset[ci_bs_key][currSM], variants_subset['dG_err'], fmax_err, fmin_err,
                                                           varyA=varyA, conc_init=conc_init, conc_init_quantile=conc_init_quantile, **kwargs)
        else:
            fitResult, predictedConc = deconvoluteMixtures(variants_subset[bs_key][currSM], variants_subset['dG'], fmax, fmin,
                                                           None, variants_subset['dG_err'], fmax_err, fmin_err,
                                                           varyA=varyA, conc_init=conc_init, conc_init_quantile=conc_init_quantile, **kwargs)
        list_predictedConc.append(predictedConc)
        dict_fitResults[currSM] = fitResult
    
    predictedConcMatrix = pd.concat(list_predictedConc, axis=1, keys=listSM).reindex(listSM)
    return predictedConcMatrix, dict_fitResults


def fitAllComplexMixtures(variants_subset, listCM, fmax=True, fmin=True, data_err=True, norm=True, 
                          varyA=False, conc_init=None, conc_init_quantile=0.1, **kwargs):
    '''Fit all complex mixtures measurements'''
    # Define column names for fmax, fmin, and bs in the normalized and unnormalzied cases
    if norm:
        fmax_key = 'fmax_norm'
        fmax_err_key = 'fmax_err_norm'
        fmin_key = 'fmin_norm'
        fmin_err_key = 'fmin_err_norm'
        cm_key = 'N'

    else:
        fmax_key = 'fmax'
        fmax_err_key = 'fmax_err'
        fmin_key = 'fmin'
        fmin_err_key = 'fmin_err'
        cm_key = 'S'
    
    # Get fmax and fmin if requested
    if fmax:
        fmax = variants_subset[fmax_key]
        fmax_err = variants_subset[fmax_err_key]
    else:
        fmax = None
        fmax_err = None
    
    if fmin:
        fmin = variants_subset[fmin_key]
        fmin_err = variants_subset[fmin_err_key]
    else:
        fmin = None
        fmin_err = None
    
    # Fit all pure samples
    list_predictedConc = []
    dict_fitResults = {}
    for currCM in listCM:
        if data_err:
            return 0
        else:
            fitResult, predictedConc = deconvoluteMixtures(variants_subset[cm_key][currCM], variants_subset['dG'], fmax, fmin,
                                                           None, variants_subset['dG_err'], fmax_err, fmin_err,
                                                           varyA=varyA, conc_init=conc_init, conc_init_quantile=conc_init_quantile, **kwargs)
        list_predictedConc.append(predictedConc)
        dict_fitResults[currCM] = fitResult
    
    predictedConcMatrix = pd.concat(list_predictedConc, axis=1, keys=listCM).reindex(listCM)
    return predictedConcMatrix, dict_fitResults


def reportFitStatusAllSamples(dict_fitResults):
    for currSM in dict_fitResults:
        print currSM+':'
        print '  ier:'+str(dict_fitResults[currSM].ier)
        print '  nfev:'+str(dict_fitResults[currSM].nfev)
        print '  lmdif_message: '+dict_fitResults[currSM].lmdif_message
        print '  message: '+dict_fitResults[currSM].message


# --- Performance metrics --- # 


def evaluatePerformance(y1, y2, metric='RMSLE', axis=0):
    if metric == 'RMSE':
        return np.sqrt(((y1 - y2)**2).mean(axis=axis))
    elif metric == 'RMSLE':
        LE = np.log(y1+1) - np.log(y2+1)
        return np.sqrt((LE**2).mean(axis=axis))
    elif metric == 'IRMSLE':
        LE = np.log(y1+1) - np.log(y2+1)
        return 1 / (1 + np.sqrt((LE**2).mean(axis=axis)))
    elif metric == 'IERMSLE':
        LE = np.log(y1+1) - np.log(y2+1)
        return np.exp(-np.sqrt((LE**2).mean(axis=axis)))
    elif metric == 'pearson':
        return y1.corrwith(y2, axis=axis)
