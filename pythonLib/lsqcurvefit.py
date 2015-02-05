# Anthony Ho, ahho@stanford.edu, 2/3/2015
# Last update 2/4/2015

import numpy as np
from scipy import optimize


class lsqcurvefit:
    """ Python class for non-linear least square fitting 
    
    This class does non-linear least squares fitting using scipy.optimize.minimize 
    
    Features of the fitting function include:
    - fitting with lower/upper bounds and arbitrary equality/inequality constrains
    - weighted fitting with errors
    - user-defined Jacobian for speed-up
    - choice of minimization algorithm (default=SLSQP). See Scipy documentation for 
      the list of methods available
    
    Arguments:
    - func needs to take the form of func(params, x)
    - params0 needs to be a list/1-D ndarray
    - x, y, sigma need to be 1-D ndarrays of the same length
    - jac needs to return a 2D-ndarray of shape=(len(params), len(x))
    - bounds needs to be a list of tuples specifying the lower and upper bound for 
      each parameter [(pl0, pu0),(pl1, pu1),...] 
    - see Scipy documentation for scipy.optimize.minimizefor the rest of the arguments

    """

    # Constructor which also does fitting
    def __init__(self, func, x, y, params0,
                 bounds=None, constraints=(), jac=None,
                 sigma=None, method='SLSQP',
                 maxiter=100, tol=None, disp=False):
        
        results = self._fit(func, x, y, params0,
                            bounds=bounds, constraints=constraints, jac=jac,
                            sigma=sigma, method=method,
                            maxiter=maxiter, tol=tol, disp=disp)
        
        self.status = results.status
        self.success = results.success
        self.residual = results.fun
        self.params = results.x
        self.message = results.message
        self.jac = results.jac
        self.nit = results.nit

    
    # Fit curves using optimize.minimize, default method = SLSQP
    def _fit(self, func, x, y, params0, bounds, constraints, jac,
            sigma, method, maxiter, tol, disp):
        
        if jac:
            jacResnorm = self._resnormPrime
        else:
            jacResnorm = None
            
        return optimize.minimize(self._resnorm, params0, args=(x, y, func, jac, sigma),
                                 bounds=bounds, constraints=(), jac=jacResnorm,
                                 method=method, 
                                 tol=tol, options={'maxiter': maxiter, 'disp': disp})

    # Compute the sum of (y_i-f(x_i))^2 or 
    # the sum of w_i*(y_i-f(x_i))^2 where w_i = 1/sigma_i^2 if sigma is given
    def _resnorm(self, params, x, y, func, funcprime, sigma=None):
        if sigma is None:
            return np.sum((func(params, x) - y)**2)
        else:
            return np.sum(((func(params, x) - y)/sigma)**2)

    # Compute the partial derivatives of the sum of (y_i-f(x_i))^2 or 
    # the partial derivatives of the sum of w_i*(y_i-f(x_i))^2 where w_i = 1/sigma_i^2 if sigma is given
    def _resnormPrime(self, params, x, y, func, funcprime, sigma=None):
        if sigma is None:
            return np.sum(2*(func(params, x) - y)*funcprime(params, x), axis=1)
        else:
            return np.sum(2*(func(params, x) - y)*funcprime(params, x)/sigma**2, axis=1)
