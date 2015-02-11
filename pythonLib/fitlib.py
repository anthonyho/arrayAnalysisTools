# Anthony Ho, ahho@stanford.edu, 2/3/2015
# Summary statistics functions adopted from Lauren Chircus's NLS class
# Last update 2/9/2015
"""Library of tools for fitting and analysis"""


import numpy as np
from scipy import optimize
from scipy import stats
import matplotlib.pyplot as plt

# To-add:
# - Log likelihood
# - Akaike criterion
# - summary function
# - plot to return axis
# - C.I.
# - bootstrapping


def fTest(model1, model2):
    if model1.DOF <= model2.DOF:
        fScore = ((model1.RSS - model2.RSS)/(float(model1.DOF) - float(model2.DOF))) / (model2.RSS/float(model2.DOF))
        pValue = stats.f.sf(fScore, model1.DOF, model2.DOF)
        return fScore, pValue
    else:
        fScore = ((model2.RSS - model1.RSS)/(float(model2.DOF) - float(model1.DOF))) / (model1.RSS/float(model1.DOF))
        pValue = stats.f.sf(fScore, model2.DOF, model1.DOF)
        return fScore, pValue


class lsqcurvefit:
    """Python class for non-linear least square fitting

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
    - jac needs to return a 2D-ndarray of shape=(len(x), len(params))
      i.e. a m x n matrix where m is the number of datapoints and n is the number of parameters
    - bounds needs to be a list of tuples specifying the lower and upper bound for
      each parameter [(pl0, pu0),(pl1, pu1),...]
    - see Scipy documentation for scipy.optimize.minimizefor the rest of the arguments
    """

    # Constructor which also does fitting
    def __init__(self, func, x, y, params0,
                 bounds=None, constraints=(), jac=None,
                 sigma=None, method='SLSQP',
                 maxiter=100, tol=None, epsilon=None, disp=True):

        # Get boolean array to indicate missing data in y
        isFiniteBoolArray = np.isfinite(y)

        # Assign instance variables
        self.func = func
        self.x = x[isFiniteBoolArray]  # Get rid of missing data in y
        self.y = y[isFiniteBoolArray]  # Get rid of missing data in y
        self.params0 = params0
        self.bounds = bounds
        self.constraints = constraints
        self.funcPrime = jac
        self.sigma = sigma
        self.method = method
        self.maxiter = maxiter
        self.tol = tol
        self.epsilon = epsilon
        self.disp = disp
        self.nDataPoints = len(y)
        self.DOF = self.nDataPoints-len(self.params0)

        # Sanity check of input parameters
        self._sanityCheck()

        # Call optimize.minimize to fit
        results = self._fit()

        # Assign some results as instance variables
        self.status = results.status
        self.success = results.success
        self.message = results.message

        self.params = results.x
        self.RSS = results.fun

        try:
            self.jacobianRSS = results.jac
        except AttributeError:
            self.jacobianRSS = None

        try:
            self.nit = results.nit
        except AttributeError:
            self.nit = None

        # Compute the Jacobian of func at params that minimizes RSS
        self.jacobianFunc = self._compute_jacobianFunc()

        # Compute fitting statistics: (see methods for formulas)

        # Compute R-squared
        self.R2 = self._compute_Rsquared()
        # Compute adjusted R-squared
        self.adjR2 = self._compute_adjustedRsquared()
        # Compute reduced Chi-squared/residual variance/mean square error
        self.reChi2 = self._compute_reducedChiSquared()
        # Compute standard error of the regression/standard error of the equation
        self.SER = self._compute_standardErrorRegression()
        # Compute standard errors of the fit parameters
        self.paramSEs = self._compute_paramSEs()
        try:
            # Compute parameters' t-statistic
            self.paramTvals = self._compute_tStatistic()
            # Compute parameters' p-values
            self.paramPvals = self._compute_pValuesFromT()
        except TypeError:
            self.paramTvals = None
            self.paramPvals = None

        # Cleanups
        del disp

    # Sanity check of input parameters
    def _sanityCheck(self):
        if not isinstance(self.x, np.ndarray) or not isinstance(self.y, np.ndarray):
            raise ValueError("x and y need to be numpy.ndarrays!")
        if len(self.x) != len(self.y):
            raise ValueError("x and y need to be the same length!")

    # Main method to fit curves using optimize.minimize, default method = SLSQP
    def _fit(self):

        if self.funcPrime:
            RSSprime = self._compute_RSSprime
        else:
            RSSprime = None

        return optimize.minimize(self._compute_RSS, self.params0,
                                 args=(self.x, self.y, self.func, self.funcPrime, self.sigma),
                                 bounds=self.bounds, constraints=self.constraints, jac=RSSprime,
                                 method=self.method,
                                 tol=self.tol, options={'maxiter': self.maxiter, 'disp': self.disp})

    # Compute the Jacobian of func at params that minimizes RSS
    # Use user-supplied Jacobian if provided, otherwise compute numerically
    def _compute_jacobianFunc(self):
        if self.funcPrime:
            return self.funcPrime(self.params, self.x)
        else:
            if self.epsilon:
                eps = self.epsilon
            else:
                eps = np.sqrt(np.finfo(np.float).eps)
            return np.array([optimize.approx_fprime(self.params, self.func, eps, xi) for xi in self.x])

    # Compute R-squared
    def _compute_Rsquared(self):
        SSres = self._compute_RSS(self.params, self.x, self.y, self.func, None, self.sigma)
        SStot = self._compute_RSS(self.params, self.x, np.mean(self.y), self.func, None, self.sigma)
        return 1 - SSres/SStot

    # Compute adjusted R-squared
    def _compute_adjustedRsquared(self):
        return 1 - (1-self.R2)*(self.nDataPoints-1)/(self.DOF)

    # Compute reduced Chi-squared/residual variance/mean square error
    #     reChiSq = sum_i[ (func(p*,x_i)-y_i)^2 ]/dof
    # in the unweighted case, or
    #     reChiSq = sum_i[ (func(p*,x_i)-y_i)^2/sigma^2 ]/dof
    # in the weighted case, where p* are the optimized parameters and dof is
    # the degrees of freedom
    def _compute_reducedChiSquared(self):
        return self.RSS/self.DOF

    # Compute the standard error of the regression/standard error of the equation
    def _compute_standardErrorRegression(self):
        return np.sqrt(self.reChi2)

    # Compute the standard error of the fitted parameters
    # The covariance matrix is given by:
    #     covar = (J^T J)^{-1}
    #
    # In the case of unweighted least square fitting, the function returns the square
    # root of the diagonal elements of reChiSq*covar, which estimates the statistical
    # error on the best-fit parameters from the scatter of the underlying data.
    #
    # In the case of weighted least square fitting, the function returns the square
    # root of the diagonal elements of covar, which estimates the statistical error on
    # the best-fit parameters resulting from the Gaussian errors sigma_i on the
    # underlying data y_i.
    #
    # Build-in except case in case of singular matrix
    def _compute_paramSEs(self):
        try:
            covar = np.linalg.inv(np.dot(self.jacobianFunc.transpose(), self.jacobianFunc))
            if self.sigma is None:
                return np.sqrt(np.diag(self.reChi2*covar))
            else:
                return np.sqrt(np.diag(covar))
        except np.linalg.linalg.LinAlgError:
            return None

    # Compute the parameter estimates' t-statistic
    def _compute_tStatistic(self):
        return self.params/self.paramSEs

    # Compute the p-values given the parameter estimates' t-statistic
    def _compute_pValuesFromT(self):
        return stats.t.sf(np.abs(self.paramTvals), self.DOF)*2

    # Objective function (residual sum of squares) to be minimized.
    # Compute sum_i[ (func(p,x_i)-y_i)^2 ] or, if sigma is given,
    # sum_i[ w_i*(func(p,x_i)-y_i)^2 ] where w_i = 1/sigma_i^2
    def _compute_RSS(self, params, x, y, func, _, sigma=None):
        if sigma is None:
            return np.sum((func(params, x) - y)**2)
        else:
            return np.sum(((func(params, x) - y)/sigma)**2)

    # Derivative of the objective function (residual sum of squares) to be minimized.
    # Compute the partial derivatives of sum_i[ (func(p,x_i)-y_i)^2 ] or, if sigma
    # is given, the partial derivatives of sum_i[ w_i*(func(p,x_i)-y_i)^2 ] where
    # w_i = 1/sigma_i^2
    def _compute_RSSprime(self, params, x, y, func, funcPrime, sigma=None):
        if sigma is None:
            return np.sum(2 * (func(params, x) - y) * funcPrime(params, x).transpose(), axis=1)
        else:
            return np.sum(2 * (func(params, x) - y) * funcPrime(params, x).transpose() / sigma**2, axis=1)

    # Make summary text in plot
    def _makeSummaryInPlot(self, paramNames=None):
        if paramNames:
            paramStrList = [u'{}: {:.4g}{}{:.4g}'.format(name, val, u'\u00B1', se) 
                            for name, val, se in zip(paramNames, self.params, self.paramSEs)]
        else:
            paramStrList = [u'param[{:d}]: {:.4g}{}{:.4g}'.format(i, val, u'\u00B1', se) 
                            for i, (val, se) in enumerate(zip(self.params, self.paramSEs))]
        return '\n'.join(paramStrList)

    # Public function to plot data and fitted curve
    def plot(self, figsize=(6, 6), numPlotPoints=500,
             markeredgecolor='r', markeredgewidth='2', markersize=10,
             linecolor='b', linewidth='2', borderwidth='2',
             xlabel=None, ylabel=None, title=None,
             summary=None, paramNames=None,
             block=False):

        # Compute the x axis points for plotting the fitted line
        xPlotPoints = np.arange(min(self.x), max(self.x)+1, (max(self.x)-min(self.x))/numPlotPoints)

        # Plot the data and fitted line
        fig = plt.figure(figsize=figsize)
        ax = plt.gca()
        fig.patch.set_facecolor('white')
        plt.plot(self.x, self.y, marker='o', linestyle='None', color='w',
                 markeredgecolor=markeredgecolor, markeredgewidth=markeredgewidth, markersize=markersize)
        plt.plot(xPlotPoints, self.func(self.params, xPlotPoints), color=linecolor, linewidth=linewidth)

        # Show labels, title and summary if requested
        if xlabel:
            ax.set_xlabel(xlabel)
        if ylabel:
            ax.set_ylabel(ylabel)
        if title:
            ax.set_title(title, y=1.02)
        if summary:
            if type(summary) != tuple:
                summary = (0.03, 0.15)
            ax.text(summary[0], summary[1], self._makeSummaryInPlot(paramNames), transform=ax.transAxes, 
                    fontsize=14, verticalalignment='top')

        # Make it pretty
        plt.rc('axes', linewidth=borderwidth)
        ax.xaxis.label.set_fontsize(16)
        ax.yaxis.label.set_fontsize(16)
        ax.title.set_fontsize(16)

        plt.show(block=block)
        return

    # Public function to print summary of the fitting statistics
    def summarize(self):
        return
