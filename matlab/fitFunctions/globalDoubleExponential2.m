function F = globalDoubleExponential2(params,X)
    T = X(:,1);        % time
    dsid = X(:,2);     % dataset id
    A0=params(1);
    A1=params(2);
    tau1 = params(3);
    tau2 = params(4:(3+max(dsid)))';
    
    
    F = A0+A1*exp(-T/tau1)+(1-A0-A1)*exp(-T./tau2(dsid));
   
    