classdef fitObj
    
    properties
        
        functionHandle
        
        params0
        paramsLB
        paramsUB
        
        constantsNames
        constantsValues
        
        options = optimset('MaxFunEvals',10000,'MaxIter',10000, ...
        'TolFun',1e-11,'display','none');
        
    end
    
    
    methods
        
        function FO = fitObj(functionHandle,params0,paramsLB,paramsUB, ...
                constantsNames,constantsValues,options)
            FO.functionHandle = functionHandle;
            FO.params0 = params0;
            FO.paramsLB = paramsLB;
            FO.paramsUB = paramsUB;
            FO.constantsNames = constantsNames;
            FO.constantsValues = constantsValues;
            if nargin >= 7
                FO.options = options;
            end
        end
        
        function [params,eof,seParams] = fitOnly(FO,xAxis,data)
            
            % Define the number of data points
            nDataPoints = length(data);
            
            % Compute the degree of freedom
            dof = nDataPoints-length(FO.params0);
            
            % Fitting using lsqcurvefit
            [params,resnorm,residual,exitflag,output,lambda,jacobian] = ...
                lsqcurvefit(FO.functionHandle,FO.params0, ...
                xAxis,data, ...
                FO.paramsLB,FO.paramsUB,FO.options);
            
            % Computing error of the fit
            eof = resnorm/dof;
            
            % Estimating error of fit parameters
            [Q,R] = qr(jacobian,0);
            mse = sum(abs(residual).^2)/(size(jacobian,1)-size(jacobian,2));
            Rinv = inv(R);
            Sigma = Rinv*Rinv'*mse;
            seParams = full(sqrt(diag(Sigma)));
            
        end
        
        function [hPoints,hLine,params,eof,seParams] = fitAndPlot(FO,xAxis,data)

            hold on;

            % Plot the data points
            hPoints = plot(xAxis,data,'o');
            
            % Calling fitOnly
            [params,eof,seParams] = FO.fitOnly(xAxis,data);
            
            % Plotting fitted line
            hLine = plot(xAxis,FO.functionHandle(params,xAxis),'-');
            
            % Make plot pretty
            makepretty;
            
        end
    
    
        
    end
    
    
    
end