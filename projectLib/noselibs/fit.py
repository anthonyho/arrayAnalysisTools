# Anthony Ho, ahho@stanford.edu, 1/5/2017
# Last update 3/9/2017
"""Python module for deconvolution of complex mixtures of ligands based on biophysical model"""


import numpy as np
import pandas as pd
import lmfit
import scipy.odr as odr
import fit_funs
import plot
import aux


class deconvoluteMixtures:
    """Python class for deconvoluting complex mixtures into individual components"""


    # Constructor which also does fitting of all samples
    def __init__(self, variants, listSamples, listLigands,
                 currConc=None, trueConcDf=None,
                 use_fmax=True, use_fmin=True,
                 use_err={'data': 1, 'dG': 1, 'fmax': 1, 'fmin': 1},
                 norm=True, varyA=False, unit='uM',
                 conc_init=None, conc_init_percentile=10,
                 method='leastsq', maxit_lmfit=0, maxit_odr=1000, **kwargs):
        """Deconvolute samples given the variant matrix"""
        # Assign instance variables
        self.listSamples = listSamples
        self.listLigands = listLigands
        self.currConc = currConc
        self.use_fmax = use_fmax
        self.use_fmin = use_fmin
        self.use_err = use_err
        self.norm = norm
        self.varyA = varyA
        self.unit = unit
        self.conc_init = conc_init
        self.conc_init_percentile = conc_init_percentile
        self.method = method
        self.maxit_lmfit = maxit_lmfit
        self.maxit_odr = maxit_odr
        self.ndata = len(variants)

        # Get fitting methods
        self._getFittingMethods()

        # Create trueConcDf if not provided
        if trueConcDf:
            self.trueConcDf = trueConcDf
        elif currConc and (listSamples == listLigands):
            self.trueConcDf = aux.generate_ligand_trueConcMat(currConc, listSamples) # to be fixed for listSamples != listLigands
        else:
            self.trueConcDf = None

        # Extract and define variables used for fitting
        self._prepareVariables(variants)

        # Fit all samples
        self.results = {}
        for sample in listSamples:
            self.results[sample] = self._fitSingleSample(sample, **kwargs)
        
        # Create predicted concentration matrix
        final_method = 'ODR' if self.run_odr else 'OLS'
        self.predictedConcMatrix = pd.concat({sample: self.results[sample][final_method].predictedConc 
                                              for sample in listSamples}, 
                                             axis=1).reindex(index=listLigands, columns=listSamples)


    # Private method to define fitting methods
    def _getFittingMethods(self):
        if self.method == 'odr':
            self.lmfit_method = None
            self.run_odr = True
        elif self.method[-4:] == '+odr':
            self.lmfit_method = self.method[:-4]
            self.run_odr = True
        else:
            self.lmfit_method = self.method
            self.run_odr = False


    # Private method to extract and define variables used for fitting
    def _prepareVariables(self, variants):
        # Define column names for fmax, fmin, data of samples, and their errors
        # in the normalized and unnormalized cases
        col = {'fmax': 'fmax',
               'fmax_err': 'fmax_err',
               'fmin': 'fmin',
               'fmin_err': 'fmin_err'}
        if self.currConc:
            col['data'] = 'bs_'+str(self.currConc)+'.0'
            col['data_err'] = 'bs_'+str(self.currConc)+'.0_err'
        else:
            col['data'] = 'cs'
            col['data_err'] = 'cs_err'
            
        if self.norm:
            col = {key: col[key]+'_norm' for key in col}

        # Extract input to be fed into the fitting algorithm
        # - data -
        self.data = variants[col['data']]
        if self.use_err['data']:
            self.data_err = variants[col['data_err']]
        else:
            self.data_err = pd.DataFrame(1, index=self.data.index, columns=self.data.columns)
        # - dG -
        self.dG = variants['dG']
        if self.use_err['dG']:
            self.dG_err = variants['dG_err']
        else:
            self.dG_err = pd.DataFrame(1, index=self.dG.index, columns=self.dG.columns)
        # - fmax -
        if self.use_fmax:
            self.fmax = variants[col['fmax']]
        else:
            self.use_err['fmax'] = 0
            self.fmax = pd.DataFrame(1, index=self.dG.index, columns=self.dG.columns)
        if self.use_err['fmax']:
            self.fmax_err = variants[col['fmax_err']]
        else: 
            self.fmax_err = pd.DataFrame(1, index=self.fmax.index, columns=self.fmax.columns)
        # - fmin -
        if self.use_fmin:
            self.fmin = variants[col['fmin']]
        else:
            self.use_err['fmin'] = 0
            self.fmin = pd.Series(0, index=self.dG.index)
        if self.use_err['fmin']:
            self.fmin_err = variants[col['fmin_err']]
        else: 
            self.fmin_err = pd.Series(1, index=self.fmin.index)
        # - ifixx -
        self.ifixx = np.ones(len(self.dG.columns) + len(self.fmax.columns) + 1)
        self.ifixx[0:len(self.dG.columns)] = self.use_err['dG']
        self.ifixx[len(self.dG.columns):len(self.dG.columns)+len(self.fmax.columns)] = self.use_err['fmax']
        self.ifixx[-1] = self.use_err['fmin']

        # Define initial concenrations for all samples
        try:
            _conc_init = np.ones(len(self.dG.columns)) * self.conc_init
            self.beta_init = aux.KdtodG(_conc_init, unit=self.unit)
        except TypeError:
            self.beta_init = np.percentile(self.dG, self.conc_init_percentile, axis=0)


    # Private method to fit a single sample
    def _fitSingleSample(self, sample, **kwargs):
        """Fit the concentrations of ligands using lmfit"""
        # Initialize and typecast variables
        x, x_err, y, y_err = \
            fit_funs._transformVariables(self.data[sample], self.dG, self.fmax, self.fmin,
                                         self.data_err[sample], self.dG_err, self.fmax_err, self.fmin_err)

        # Construct parameters
        params = lmfit.Parameters()
        params.add('A', value=1.0, min=0.0, vary=self.varyA)
        for i, ligand in enumerate(self.dG.columns):
            params.add(ligand, value=self.beta_init[i])

        # Fit and extract parameter estimates
        result = {}
        if self.lmfit_method is not None:
            lmfit_result = lmfit.minimize(fit_funs._switchingEq_residuals_lmfit, params,
                                          args=(x, x_err, y, y_err, self.use_err),
                                          method=self.lmfit_method, maxfev=self.maxit_lmfit, **kwargs)
            lmfit_result.predictedConc = aux.dGtoKd(pd.Series(lmfit_result.params.valuesdict()).drop('A'), 
                                                    unit=self.unit)
            result['OLS'] = lmfit_result
        if self.run_odr:
            try:
                beta0 = np.array(pd.Series(lmfit_result.params.valuesdict()).drop('A'))
            except NameError:
                beta0 = self.beta_init
            odr_model = odr.Model(fit_funs._switchingEq,
                                  fjacb=fit_funs._switchingEq_jacobian_mu,
                                  fjacd=fit_funs._switchingEq_jacobian_x)
            odr_data = odr.RealData(x, y, sx=x_err, sy=y_err)
            odr_obj = odr.ODR(odr_data, odr_model, beta0=beta0, ifixx=self.ifixx, maxit=self.maxit_odr)
            odr_result = odr_obj.run()
            odr_result.predictedConc = aux.dGtoKd(pd.Series(odr_result.beta, index=self.dG.columns), unit='uM')
            result['ODR'] = odr_result

        return result

















    # Public method to plot the predicted concentration confusion matrix
    def plotPredictedConcMatrix(self, setup='', fig_dir=None):
        """Plot the predicted concentration confusion matrix"""
        cg = plot.plotPredictedConcMatrix(fitResults=self, setup=setup, fig_dir=fig_dir)
        return cg


    # Public method to plot the convergence and performance metrics
    def plotFitStatus(self, setup='', metric='IERMSLE', fig_dir=None):
        """Plot the convergence and performance metrics"""
        plot.plotFitStatus(fitResults=self, setup=setup, metric=metric, fig_dir=fig_dir)


    # Public method to report the fit status of all samples
    def reportFitStatus(self):
        """Report the status of the fits"""
        for sample in self.results:
            trueRedChi = self.reportFit(sample, output='redchi', weighted=True, params='true')
            print '{:<6} ier = {}, nfev = {}, redchi = {:.3f}, trueRedChi = {:.3f}'.format(sample+':', 
                                                                                           self.results[sample].ier, 
                                                                                           self.results[sample].nfev, 
                                                                                           self.results[sample].redchi, 
                                                                                           trueRedChi)
            print '       message: {}'.format(self.results[sample].message.replace('\n', ''))
            if self.method == 'leastsq':
                print '       lmdif_message: {}'.format(self.results[sample].lmdif_message.replace('\n', ''))


    # Public method to evaluate the performance of the fit
    def evaluatePerformance(self, metric='IERMSLE', axis=0):
        """Evaulate the performance of the fits"""
        # Unpack variables
        y1 = self.predictedConcMatrix
        y2 = self.trueConcDf

        # Evaluate performance
        if metric == 'RMSE':
            RMSE = np.sqrt(((y1 - y2)**2).mean(axis=axis))
            return RMSE.dropna()
        elif metric == 'RMSLE':
            LE = np.log(y1+1) - np.log(y2+1)
            RMSLE = np.sqrt((LE**2).mean(axis=axis))
            return RMSLE.dropna()
        elif metric == 'IRMSLE':
            LE = np.log(y1+1) - np.log(y2+1)
            IRMSLE = 1 / (1 + np.sqrt((LE**2).mean(axis=axis)))
            return IRMSLE.dropna()
        elif metric == 'IERMSLE':
            LE = np.log(y1+1) - np.log(y2+1)
            IERMSLE = np.exp(-np.sqrt((LE**2).mean(axis=axis)))
            return IERMSLE.dropna()
        elif metric == 'pearson':
            pearson = y1.corrwith(y2, axis=axis)
            return pearson.dropna()


    # Public method to compute various functional terms given parameters
    def reportFit(self, sample, output, weighted=True, params='fitted'):
        """Return various functional terms using the fitted parameters"""
        # Initialize and typecast variables
        _data, _dG, _fmax, _fmin, \
        _data_err, _dG_err, _fmax_err, _fmin_err = \
        fit_funs._prepareVariables(self.data[sample], self.dG, self.fmax, self.fmin,
                                   self.data_err[sample], self.dG_err, self.fmax_err, self.fmin_err)

        # Extract concentrations of ligands from params dict
        if params == 'fitted':
            A, mu = fit_funs._extractParams(self.results[sample].params)
        elif params == 'true':
            A = 1
            mu = [aux.KdtodG(self.trueConcDf[sample][ligand], unit=self.unit) 
                  for ligand in self.dG.columns]
            mu = np.array(mu).reshape(1, -1)
        elif params == 'init':
            A = 1
            mu = self.mu_init
        elif params == 'zeros':
            A = 1
            _conc = np.zeros(len(self.dG.columns))
            mu = aux.KdtodG(_conc, unit=self.unit).reshape(1, -1)
        else:
            try:
                A = 1
                _conc = np.ones(len(self.dG.columns)) * params
                mu = aux.KdtodG(_conc, unit=self.unit).reshape(1, -1)
            except TypeError:
                A, mu = fit_funs._extractParams(params)

        # Define weights if requested
        if weighted:
            weights = fit_funs._switchingEq_errors(mu, _dG, _fmax, _fmin, _data_err, _dG_err, _fmax_err, _fmin_err, self.use_err)
        else:
            weights = 1

        # Output
        if output=='data':
            return _data / weights
        elif output=='eq':
            return A * fit_funs._switchingEq(mu, _dG, _fmax, _fmin) / weights
        elif output=='residual':
            return (_data - A * fit_funs._switchingEq(mu, _dG, _fmax, _fmin)) / weights
        elif output=='chisqr':
            r = (_data - A * fit_funs._switchingEq(mu, _dG, _fmax, _fmin)) / weights
            return np.linalg.norm(r)**2
        elif output=='redchi':
            r = (_data - A * fit_funs._switchingEq(mu, _dG, _fmax, _fmin)) / weights
            v = _dG.shape[0] - _dG.shape[1]
            return np.linalg.norm(r)**2 / v
        elif output=='err':
            return fit_funs._switchingEq_errors(mu, _dG, _fmax, _fmin, _data_err, _dG_err, _fmax_err, _fmin_err, self.use_err)
        elif output=='err_mu':
            return fit_funs._switchingEq_jacobian_mu(mu, _dG, _fmax, _fmin)
        elif output=='err_dG':
            return fit_funs._switchingEq_jacobian_dG(mu, _dG, _fmax, _fmin)
        elif output=='err_fmax':
            return fit_funs._switchingEq_jacobian_fmax(mu, _dG)
        elif output=='err_fmin':
            return fit_funs._switchingEq_jacobian_fmin(mu, _dG)
