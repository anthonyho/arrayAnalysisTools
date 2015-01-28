function allFitParams = singleClusterFit(data,xAxis,functionHandle, ...
    fitParams0,fitParamsLB,fitParamsUB,makefig)
% SINGLECLUSTERFIT Fits single clusters
%
% Anthony Ho, 10/27/2014


    %% Constants and initializations

    allFitParams = [];
    nClusters = size(data,1);
    
    options = optimset('MaxFunEvals',10000,'MaxIter',10000, ...
        'TolFun',1e-11,'display','none');
    
    dof = size(data,2)-length(fitParams0);
    
    
    %% Fitting
    
    for i = 1:nClusters
        
        if mod(i,100)==0
            disp(sprintf('      Fitting the %dth cluster...',i));
        end
        
        [fitParams,resnorm,residual,exitflag,output,lambda,jacobian] ...
            = lsqcurvefit(functionHandle,fitParams0, ...
            xAxis,data(i,:)',fitParamsLB,fitParamsUB,options);
        
        % Error of fit
        fitParams = [fitParams; resnorm/dof];
        
        % Estimating error of fit parameters
        [Q,R] = qr(jacobian,0);
        mse = sum(abs(residual).^2)/(size(jacobian,1)-size(jacobian,2));
        Rinv = inv(R);
        Sigma = Rinv*Rinv'*mse;
        se = full(sqrt(diag(Sigma)));

        fitParams = [fitParams; se];
        
        % Appending to the list of fit parameters
        allFitParams = [allFitParams; fitParams'];
        
        if makefig
        
            h1 = plot(xAxis,data(i,:)','ob');
            h2 = plot(xAxis,functionHandle(fitParams,xAxis),'-b');
            
            % Add display for fit parameters
            
            makepretty;
            
            k = waitforbuttonpress;
            
            delete(h1);
            delete(h2);
            
        end

    end
    
    
end