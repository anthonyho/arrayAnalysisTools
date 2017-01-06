# Anthony Ho, ahho@stanford.edu, 1/5/2016
# Last update 1/5/2015
"""Library containing the switching equations, their derivatives, residuals, and fitting functions"""


import numpy as np


def switchingEq(conc, Kd, fmax, fmin):
    '''Compute the intensity of a switching aptamer'''
    R = conc / Kd
    R_sum = R.sum(axis=1).reshape(-1, 1)
    R_weighted_sum = (R * (fmax - fmin)).sum(axis=1).reshape(-1, 1)
    return (R_weighted_sum / (1 + R_sum) + fmin).reshape(-1)


def switchingEq_jaconbian_conc(conc, Kd, fmax, fmin):
    '''Compute the Jacobian of the switching equation with respect to conc'''
    R = conc / Kd
    R_sum = R.sum(axis=1).reshape(-1, 1)
    R_weighted_sum = (R * (fmax - fmin)).sum(axis=1).reshape(-1, 1)
    return ((fmax - fmin) / Kd * (1 + R_sum) - 1 / Kd * R_weighted_sum) / (1 + R_sum)**2


def switchingEq_jacobian_Kd(conc, Kd, fmax, fmin):
    '''Compute the Jacobian of the switching equation with respect to Kd'''
    R = conc / Kd
    R_sum = R.sum(axis=1).reshape(-1, 1)
    R_weighted_sum = (R * (fmax - fmin)).sum(axis=1).reshape(-1, 1)
    return (-(fmax - fmin) * conc / (Kd**2) * (1 + R_sum) + conc / (Kd**2) * R_weighted_sum) / (1 + R_sum)**2


def switchingEq_jacobian_fmax(conc, Kd):
    '''Compute the Jacobian of the switching equation with respect to fmax'''
    R = conc / Kd
    R_sum = R.sum(axis=1).reshape(-1, 1)
    return R / (1 + R_sum)


def switchingEq_jacobian_fmin(conc, Kd):
    '''Compute the Jacobian of the switching equation with respect to fmin'''
    R = conc / Kd
    R_sum = R.sum(axis=1).reshape(-1, 1)
    return 1 / (1 + R_sum)


def switchingEq_errors(conc, Kd, fmax, fmin, Kd_err, fmax_err, fmin_err):
    '''Compute the total uncertainty of each aptamer'''
    var_Kd_sum = ((switchingEq_jacobian_Kd(conc, Kd, fmax, fmin) * Kd_err)**2).sum(axis=1)
    var_fmax_sum = ((switchingEq_jacobian_fmax(conc, Kd) * fmax_err)**2).sum(axis=1)
    var_fmin = ((switchingEq_jacobian_fmin(conc, Kd) * fmin_err)**2).reshape(-1)
    return np.sqrt(var_Kd_sum + var_fmax_sum + var_fmin)


def switchingEq_residuals(params, data, Kd, fmax, fmin, 
                          Kd_err=None, fmax_err=None, fmin_err=None):
    '''Compute the residuals of each aptamer'''
    # Extract concentrations of ligands from params dict
    if isinstance(params, lmfit.Parameters):
        conc = np.zeros((1, len(Kd.columns)))
        paramvals = params.valuesdict()
        for i, currSM in enumerate(Kd.columns):
            conc[0, i] = paramvals[currSM]
    else:
        conc = np.array(params).reshape(1, -1)
    
    # Typecast and reshape input variables
    _data = np.array(data)
    _Kd = np.array(Kd)
    _fmax = np.array(fmax)
    _fmin = np.array(fmin).reshape(-1, 1)
    _Kd_err = np.array(Kd_err)
    _fmax_err = np.array(fmax_err)
    _fmin_err = np.array(fmin_err).reshape(-1, 1)
    
    # Reture residuals of each aptamers 
    if Kd_err is None:
        return _data - paramvals['A'] * switchingEq(conc, _Kd, _fmax, _fmin)
    else:
        return (_data - paramvals['A'] * switchingEq(conc, _Kd, _fmax, _fmin)) / switchingEq_errors(conc, _Kd, _fmax, _fmin, _Kd_err, _fmax_err, _fmin_err)


def deconvoluteMixtures(I, Kd, fmax=None, fmin=None, 
                        Kd_err=None, fmax_err=None, fmin_err=None, 
                        varyA=False, conc_init=None, **kwargs):
    '''Fit the concentrations of ligands using lmfit'''
    # Define fmax and fmin and their respective errors if not provided
    if fmax is None:
        fmax = np.ones(Kd.shape)
    if fmin is None:
        fmin = np.zeros(len(Kd))
    if fmax_err is None:
        fmax_err = np.zeros(Kd.shape)
    if fmin_err is None:
        fmin_err = np.zeros(len(Kd))
       
    # Define concenrations
    if conc_init is None:
        conc_init = np.ones(len(Kd.columns)) * 0.1
    elif isinstance(conc_init, int) or isinstance(conc_init, float):
        conc_init = np.ones(len(Kd.columns)) * conc_init
    params = lmfit.Parameters()
    params.add('A', value=1.0, min=0.0, vary=varyA)
    for i, currSM in enumerate(Kd.columns):
        params.add(currSM, value=conc_init[i], min=0.0)
    
    # Fit and extract params
    fitResults = lmfit.minimize(switchingEq_residuals, params,
                                args=(I, Kd, fmax, fmin, Kd_err, fmax_err, fmin_err), 
                                **kwargs)
    predictedConc = pd.Series(fitResults.params.valuesdict()).drop('A')
    
    return fitResults, predictedConc


def fitAllPureSamples(variants_subset, currConc, listSM, fmax=True, fmin=True, norm=True, 
                      varyA=False, conc_init=None, **kwargs):
    '''Fit all the pure sample measurements at a given concentration'''
    # Define column names for fmax, fmin, and bs in the normalized and unnormalzied cases
    if norm:
        fmax_key = 'fmax_norm'
        fmax_err_key = 'fmax_err_norm'
        fmin_key = 'fmin_norm'
        fmin_err_key = 'fmin_err_norm'
        bs_key = 'bs_'+str(currConc)+'.0_norm'
    else:
        fmax_key = 'fmax'
        fmax_err_key = 'fmax_err'
        fmin_key = 'fmin'
        fmin_err_key = 'fmin_err'
        bs_key = 'bs_'+str(currConc)+'.0'
    
    # Get fmax and fmin if requested
    if fmax:
        fmax = variants_subset[fmax_key]
        fmax_err = variants_subset[fmax_err_key]
    else:
        fmax = None
    
    if fmin:
        fmin = variants_subset[fmin_key]
        fmin_err = variants_subset[fmin_err_key]
    else:
        fmin = None
        
    # Fit all pure samples
    list_predictedConc = []
    dict_fitResults = {}
    for currSM in listSM:
        fitResults, predictedConc = deconvoluteMixtures(variants_subset[bs_key][currSM], variants_subset['Kd'], fmax, fmin,
                                                        variants_subset['Kd_err'], fmax_err, fmin_err,  
                                                        varyA=varyA, conc_init=None, **kwargs)
        list_predictedConc.append(predictedConc)
        dict_fitResults[currSM] = fitResults
    
    predictedConcMatrix = pd.concat(list_predictedConc, axis=1, keys=listSM).reindex(listSM)
    return predictedConcMatrix, dict_fitResults
