# Anthony Ho, ahho@stanford.edu, 2/3/2015
# Summary statistics functions adopted from Lauren Chircus's NLS class
# Last update 2/11/2015
"""Library of tools for fitting and analysis"""


import numpy as np
from scipy import optimize
from scipy import stats
import matplotlib.pyplot as plt

# To-add:
# - Log likelihood
# - Akaike criterion
# - C.I.
# - bootstrapping


class lsqcurvefit:
    """Python class for non-linear least square fitting using scipy.optimize.minimize

       Features of this fitting class include:
       - fitting with lower/upper bounds and arbitrary equality/inequality constrains
       - weighted least square fitting (weight = 1/sigma^2)
       - user-defined Jacobian for speed-up
       - choice of minimization algorithm (default=SLSQP). See Scipy documentation for the list
         of methods available
       - automatically handle missing data (nan, inf, -inf)
       - auxiliary functions to show fitted curve against datapoints and print summary of 
         fitting statistics

    Usages:

    To fit:

       fitObj = fitlib.lsqcurvefit(func, x, y, params0,
                                   bounds=None, constraints=(), jac=None,
                                   sigma=None, method='SLSQP',
                                   maxiter=200, tol=None, epsilon=None, disp=True)

       Required arguments:
       - func: callable
           Model function to be fitted to the data. Must take the form of func(params, x)
       - x: M-length sequence
           Independent variable to be passed to func()
       - y: M-length sequence
           Dependent data 
       - params0: N-length sequence
           Initial guess for the parameters

       Optional arguments:
       - bounds: N-length sequence of tuples
           Bounds for variables (only for L-BFGS-B, TNC and SLSQP). A list of tuples specifying 
           the lower and upper bound for each parameter [(pl0, pu0),(pl1, pu1),...]. Use None 
           for one of min or max when there is no bound in that direction
           (default=None)
       - constraints: dict or sequence of dict
           Constraints definition (only for COBYLA and SLSQP). See documentation for 
           scipy.optimize.minimize for syntax
           (default=())
       - jac: callable
           Jacobian (gradient) of model function. Must take the same arguments as func() and 
           return a 2D-ndarray of shape=(M, N) where M is the number of datapoints and N is the 
           number of parameters
           (default=None)
       - sigma: M-length sequence
           If provided, these values are used as weights in the least-squares problem as 
           weight = 1/sigma^2
           (default=None)
       - method: str
           Type of solver. See documentation for scipy.optimize.minimize for details
           (default='SLSQP')
       - maxiter: int
           Maximum number of iterations to perform
           (default=200)
       - tol: float
           Tolerance for termination
       - epsilon: float
           Step size used for numerical approximation of the jacobian of the model function
       - disp: bool
           Set to True to print convergence messages
           (default=True)

    To plot:

       ax = fitObj.plot(figsize=(7.5, 7.5), numPlotPoints=500,
                        markeredgecolor='r', markeredgewidth='2', markersize=10,
                        linecolor='b', linewidth='2', borderwidth='2',
                        xlabel=None, ylabel=None, title=None,
                        summary=None, paramNames=None, block=False)

       Optional arguments:
       - figsize: (w,h) tuple in inches
       - numPlotPoints: number of points to use for plotting fitted curve
       - markeredgecolor: color of the datapoints
       - markeredgewidth: edge width of the datapoints
       - markersize: size of the datapoints
       - linecolor: color of the fitted curve
       - linewidth: width of the fitted curve
       - borderwidth: width of the border in the plot
       - xlabel: x label
       - ylaebl: y label
       - title: title
       - summary: bool or tuple. If true or tuple, show fitting statistics summary in plot
                  If tuple, show the summary at the indicated location (in relative coord.)
       - paramNames: names of the parameters to display 
       - block: plotting without blocking execution

    To print summary:

       fitObj.printSummary(paramNames=None)

       Optional arguments:
       - paramNames: N-length sequence of str
           Names of the parameters to display
    """

    # Constructor which also does fitting
    def __init__(self, func, x, y, params0,
                 bounds=None, constraints=(), jac=None,
                 sigma=None, method='SLSQP',
                 maxiter=200, tol=None, epsilon=None, disp=True):

        # Get boolean array to indicate missing data in y
        x_np = np.array(x)
        y_np = np.array(y)
        isFiniteBoolArray = np.isfinite(y_np)

        # Assign instance variables
        self.func = func
        self.x = x_np[isFiniteBoolArray]  # Get rid of missing data in y
        self.y = y_np[isFiniteBoolArray]  # Get rid of missing data in y
        self.params0 = np.array(params0)
        self.bounds = bounds
        self.constraints = constraints
        self.funcPrime = jac
        if sigma:
            self.sigma = np.array(sigma)
        else: 
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

        del disp

    # Sanity check of input parameters
    def _sanityCheck(self):
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

    # Make summary text in plot
    def _makeSummaryInPlot(self, paramNames=None):
        if paramNames:
            paramStrList = [u"{}: {:.4g} {} {:.4g}".format(name, val, u'\u00B1', se)
                            for name, val, se in zip(paramNames, self.params, self.paramSEs)]
        else:
            paramStrList = [u"p[{:d}]: {:.4g} {} {:.4g}".format(i, val, u'\u00B1', se)
                            for i, (val, se) in enumerate(zip(self.params, self.paramSEs))]
        paramStrList.append(r"$\chi^2_{red}$: " + "{: .4g}".format(self.reChi2))
        paramStrList.append("SER: " + "{: .4g}".format(self.SER))
        return '\n'.join(paramStrList)

    # Public function to plot data and fitted curve
    def plot(self, figsize=(7.5, 7.5), numPlotPoints=500,
             markeredgecolor='r', markeredgewidth='2', markersize=10,
             linecolor='b', linewidth='2', borderwidth='2',
             xlabel=None, ylabel=None, title=None,
             summary=None, paramNames=None, block=False):
        """Plot the fitted curve against the datapoints """
        # Compute the x axis points for plotting the fitted line
        xPlotPoints = np.arange(min(self.x), max(self.x)+1, (max(self.x)-min(self.x))/numPlotPoints)

        # Plot the data and fitted line
        fig = plt.figure(figsize=figsize)
        fig.patch.set_facecolor('w')
        plt.plot(self.x, self.y, marker='o', linestyle='None', color='w',
                 markeredgecolor=markeredgecolor, markeredgewidth=markeredgewidth, markersize=markersize)
        plt.plot(xPlotPoints, self.func(self.params, xPlotPoints), color=linecolor, linewidth=linewidth)
        ax = plt.gca()

        # Show labels, title and summary if requested
        if xlabel:
            ax.set_xlabel(xlabel)
        if ylabel:
            ax.set_ylabel(ylabel)
        if title:
            ax.set_title(title, y=1.02)
        if summary:
            if type(summary) != tuple:
                summary = (0.98, 0.98)
            ax.text(summary[0], summary[1], self._makeSummaryInPlot(paramNames), transform=ax.transAxes,
                    fontsize=14, verticalalignment='top', horizontalalignment='right')
            self.printSummary(paramNames)

        # Make it pretty
        plt.rc('axes', linewidth=borderwidth)
        ax.xaxis.label.set_fontsize(16)
        ax.yaxis.label.set_fontsize(16)
        ax.title.set_fontsize(16)

        plt.show(block=block)

        return ax

    # Public function to print summary of the fitting statistics
    def printSummary(self, paramNames=None):
        """Print a summary of the fitting statistics"""
        # Make parameter names if not provided
        if not paramNames:
            paramNames = ["p[{:d}]".format(i) for i in range(len(self.params))]

        # Define row format for tabular output
        rowFormatHeader = "{:>15}" * 6
        rowFormatBody = "{:>15}" + "{:>15.6g}" * 4 + "{:>15}"

        # Showing summary
        print "Fitted data with {}\n".format(self.func.__name__)

        print rowFormatHeader.format("Parameter", "Estimate", "Std. error",
                                     "t-statistic", "p-value", "bounds")
        for name, est, se, tval, pval, bound in zip(paramNames, self.params,
                                                    self.paramSEs, self.paramTvals,
                                                    self.paramPvals, self.bounds):
            print rowFormatBody.format(name, est, se, tval, pval, bound)

        print "\nDegrees of freedom = {:d}".format(self.DOF)
        print "Residual sum of squares = {:.6g}".format(self.RSS)
        print "Reduced Chi-squared/residual variance/mean square error = {:.6g}".format(self.reChi2)
        print "Standard error of regression = {:.6g}".format(self.SER)
        print "R-squared = {:.6g}".format(self.R2)
        print "Adjusted R-squared = {:.6g}".format(self.adjR2)
        print "Number of iterations to convergence = {:d}".format(self.nit)

        return


# This function does NOT belong to the lsqcurvefit class
def fTest(model1, model2):
    """Do a F-test against two lsqcurvefit objects fitted with different models and return F score and p value"""
    if model1.DOF <= model2.DOF:
        fScore = ((model1.RSS - model2.RSS)/(float(model1.DOF) - float(model2.DOF))) / (model2.RSS/float(model2.DOF))
        pValue = stats.f.sf(fScore, model1.DOF, model2.DOF)
        return fScore, pValue
    else:
        fScore = ((model2.RSS - model1.RSS)/(float(model2.DOF) - float(model1.DOF))) / (model1.RSS/float(model1.DOF))
        pValue = stats.f.sf(fScore, model2.DOF, model1.DOF)
        return fScore, pValue
