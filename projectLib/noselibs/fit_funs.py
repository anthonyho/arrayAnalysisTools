# Anthony Ho, ahho@stanford.edu, 1/5/2017
# Last update 3/13/2017
"""Library containing the switching equation, its derivatives, and residuals"""


import numpy as np
import lmfit
import aux


RT = aux.RT


def _switchingEq(beta, x):
    """Compute the intensity of a switching aptamer. Shape=(m, n)"""
    mu, dG, fmax, fmin = _unpack_variables(beta, x)
    B = np.exp(-(dG - mu)/RT)
    Q = B.sum(axis=0)
    Q_weighted = (B * (fmax - fmin)).sum(axis=0)
    return (Q_weighted / (1 + Q) + fmin)


def _switchingEq_jacobian_mu(beta, x):
    """Compute the Jacobian of the switching equation
    with respect to mu. Shape=(m, n)"""
    mu, dG, fmax, fmin = _unpack_variables(beta, x)
    B = np.exp(-(dG - mu)/RT)
    Q = B.sum(axis=0)
    Q_weighted = (B * (fmax - fmin)).sum(axis=0)
    return B / RT * ((fmax - fmin) * (1 + Q) - Q_weighted) / (1 + Q)**2


def _switchingEq_jacobian_dG(beta, x):
    """Compute the Jacobian of the switching equation
    with respect to dG. Shape=(m, n)"""
    mu, dG, fmax, fmin = _unpack_variables(beta, x)
    B = np.exp(-(dG - mu)/RT)
    Q = B.sum(axis=0)
    Q_weighted = (B * (fmax - fmin)).sum(axis=0)
    return - B / RT * ((fmax - fmin) * (1 + Q) - Q_weighted) / (1 + Q)**2


def _switchingEq_jacobian_fmax(beta, x):
    """Compute the Jacobian of the switching equation
    with respect to fmax. Shape=(m, n)"""
    mu, dG, fmax, fmin = _unpack_variables(beta, x)
    B = np.exp(-(dG - mu)/RT)
    Q = B.sum(axis=0)
    return B / (1 + Q)


def _switchingEq_jacobian_fmin(beta, x):
    """Compute the Jacobian of the switching equation
    with respect to fmin. Shape=(n,)"""
    mu, dG, fmax, fmin = _unpack_variables(beta, x)
    B = np.exp(-(dG - mu)/RT)
    Q = B.sum(axis=0)
    return 1 / (1 + Q)


def _switchingEq_jacobian_x(beta, x):
    """Compute the Jacobian of the switching equation
    with respect to dG, fmax, and fmin. Shape=(2m+1, n)"""
    return np.vstack([_switchingEq_jacobian_dG(beta, x),
                      _switchingEq_jacobian_fmax(beta, x),
                      _switchingEq_jacobian_fmin(beta, x)])


def _switchingEq_errors(beta, x, x_err, y_err, use_err):
    """Compute the total uncertainty of each aptamer. Shape=(n,)"""
    dG_err, fmax_err, fmin_err, data_err = _unpack_uncertainties(x_err, y_err)
    var_data = data_err**2
    var_dG_sum = ((_switchingEq_jacobian_dG(beta, x) * dG_err)**2).sum(axis=0)
    var_fmax_sum = ((_switchingEq_jacobian_fmax(beta, x) * fmax_err)**2).sum(axis=0)
    var_fmin = ((_switchingEq_jacobian_fmin(beta, x) * fmin_err)**2)
    if np.sum(use_err) == 0:
        return 1
    else:
        return np.sqrt(var_data * use_err['data'] +
                       var_dG_sum * use_err['dG'] +
                       var_fmax_sum * use_err['fmax'] +
                       var_fmin * use_err['fmin'])


def _switchingEq_residuals(beta, x, x_err, y, y_err, A, use_err):
    """Compute the residuals of each aptamer"""
    return (y - A * _switchingEq(beta, x)) / _switchingEq_errors(beta, x, x_err, y_err, use_err)


def _switchingEq_residuals_lmfit(params, x, x_err, y, y_err, use_err):
    """Compute the residuals of each aptamer"""
    A, beta = _extract_params(params)
    return _switchingEq_residuals(beta, x, x_err, y, y_err, A, use_err)


def _unpack_variables(beta, x, y=None):
    """Unpack independent variable x into dG, fmax, and fmin"""
    mu = beta.reshape(-1, 1)
    nLigands = (len(x) - 1)/2
    dG = x[0:nLigands]
    fmax = x[nLigands:-1]
    fmin = x[-1]
    if y is None:
        return mu, dG, fmax, fmin
    else:
        data = y
        return mu, dG, fmax, fmin, data


def _unpack_uncertainties(x_err, y_err=None):
    """Unpack errors on independent variable x
    into dG_err, fmax_err, and fmin_err"""
    nLigands = (len(x_err) - 1)/2
    dG_err = x_err[0:nLigands]
    fmax_err = x_err[nLigands:-1]
    fmin_err = x_err[-1]
    if y_err is None:
        return dG_err, fmax_err, fmin_err
    else:
        data_err = y_err
        return dG_err, fmax_err, fmin_err, data_err


def _extract_params(params):
    """Extract concentrations of ligands from params dict
    into A and beta (m,)"""
    if isinstance(params, lmfit.Parameters):
        paramvals = params.valuesdict()
        A = paramvals.pop('A')
        beta = np.zeros(len(paramvals))
        for i, ligand in enumerate(paramvals):
            beta[i] = paramvals[ligand]
    else:
        A = params[0]
        beta = np.array(params[1:]).reshape(-1)
    return A, beta


def _transform_variables(data, dG, fmax, fmin,
                         data_err, dG_err, fmax_err, fmin_err):
    """Initialize and typecast variables"""
    # Typecast and reshape input variables
    _data = np.array(data)
    _dG = np.array(dG)
    _fmax = np.array(fmax)
    _fmin = np.array(fmin).reshape(-1, 1)
    _data_err = np.array(data_err)
    _dG_err = np.array(dG_err)
    _fmax_err = np.array(fmax_err)
    _fmin_err = np.array(fmin_err).reshape(-1, 1)

    # Pack variables
    x = np.hstack([_dG, _fmax, _fmin]).transpose()
    x_err = np.hstack([_dG_err, _fmax_err, _fmin_err]).transpose()
    y = _data
    y_err = _data_err

    return x, x_err, y, y_err
