function F = globalDoubleExponential3(params,X)
    T = X(:,1);        % time
    dsid = X(:,2);     % dataset id
    A0=params(1);
    tau1 = params(2);
    A2 = params(3:(2+max(dsid)))';
    tau2 = params((3+max(dsid)):(2+2*max(dsid)))';
    
    
    F = A0+(1-A0-A2(dsid)).*exp(-T/tau1)+A2(dsid).*exp(-T./tau2(dsid));
   
    