# Anthony Ho, ahho@stanford.edu, 1/5/2017
# Last update 2/7/2017
"""Library containing the switching equation, its derivatives, and residuals"""


import numpy as np
import lmfit
import aux


RT = aux.RT


def _switchingEq(mu, dG, fmax, fmin):
    """Compute the intensity of a switching aptamer"""
    B = np.exp(-(dG - mu)/RT)
    Q = B.sum(axis=1).reshape(-1, 1)
    Q_weighted = (B * (fmax - fmin)).sum(axis=1).reshape(-1, 1)
    return (Q_weighted / (1 + Q) + fmin).reshape(-1)


def _switchingEq_jacobian_mu(mu, dG, fmax, fmin):
    """Compute the Jacobian of the switching equation with respect to conc"""
    B = np.exp(-(dG - mu)/RT)
    Q = B.sum(axis=1).reshape(-1, 1)
    Q_weighted = (B * (fmax - fmin)).sum(axis=1).reshape(-1, 1)
    return B / RT * ((fmax - fmin) * (1 + Q) - Q_weighted) / (1 + Q)**2


def _switchingEq_jacobian_dG(mu, dG, fmax, fmin):
    """Compute the Jacobian of the switching equation with respect to dG"""
    B = np.exp(-(dG - mu)/RT)
    Q = B.sum(axis=1).reshape(-1, 1)
    Q_weighted = (B * (fmax - fmin)).sum(axis=1).reshape(-1, 1)
    return - B / RT * ((fmax - fmin) * (1 + Q) - Q_weighted) / (1 + Q)**2


def _switchingEq_jacobian_fmax(mu, dG):
    """Compute the Jacobian of the switching equation with respect to fmax"""
    B = np.exp(-(dG - mu)/RT)
    Q = B.sum(axis=1).reshape(-1, 1)
    return B / (1 + Q)


def _switchingEq_jacobian_fmin(mu, dG):
    """Compute the Jacobian of the switching equation with respect to fmin"""
    B = np.exp(-(dG - mu)/RT)
    Q = B.sum(axis=1).reshape(-1, 1)
    return 1 / (1 + Q)


def _switchingEq_errors(mu, dG, fmax, fmin, data_err, dG_err, fmax_err, fmin_err, use_err):
    """Compute the total uncertainty of each aptamer"""
    var_data = (data_err**2).reshape(-1)
    var_dG_sum = ((_switchingEq_jacobian_dG(mu, dG, fmax, fmin) * dG_err)**2).sum(axis=1)
    var_fmax_sum = ((_switchingEq_jacobian_fmax(mu, dG) * fmax_err)**2).sum(axis=1)
    var_fmin = ((_switchingEq_jacobian_fmin(mu, dG) * fmin_err)**2).reshape(-1)
    if np.sum(use_err) == 0:
        return 1
    else:
        return np.sqrt(var_data * use_err[0] + var_dG_sum * use_err[1] + 
                       var_fmax_sum * use_err[2] + var_fmin * use_err[3])


def _switchingEq_residuals(params, data, dG, fmax, fmin,
                           data_err, dG_err=None, fmax_err=None, fmin_err=None, use_err=[1, 1, 1, 1]):
    """Compute the residuals of each aptamer"""
    # Extract concentrations of ligands from params dict
    A, mu = _extractParams(params)

    # Reture residuals of each aptamers 
    if dG_err is None:
        return data - A * _switchingEq(mu, dG, fmax, fmin)
    else:
        return (data - A * _switchingEq(mu, dG, fmax, fmin)) / _switchingEq_errors(mu, dG, fmax, fmin, data_err, dG_err, fmax_err, fmin_err, use_err)


def _extractParams(params):
    """Extract concentrations of ligands from params dict"""
    if isinstance(params, lmfit.Parameters):
        paramvals = params.valuesdict()
        A = paramvals.pop('A')
        mu = np.zeros((1, len(paramvals)))
        for i, ligand in enumerate(paramvals):
            mu[0, i] = paramvals[ligand]
    else:
        A = params[0]
        mu = np.array(params[1:]).reshape(1, -1)
    return A, mu


def _prepareVariables(data, dG, fmax=None, fmin=None,
                      data_err=None, dG_err=None, fmax_err=None, fmin_err=None):
    """Initialize and typecast variables"""
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
