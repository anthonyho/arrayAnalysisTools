function F = globalDoubleExponential(params,X)
    T = X(:,1);        % time
    dsid = X(:,2);     % dataset id
    A0 = params(1:max(dsid))';          % different A0 for each datasets
    A1 = params(1+max(dsid));           % same A1 for all datasets
    tau1 = params(2+max(dsid));         % same tau1 for all datasets
    A2 = params((3+max(dsid)):(2+2*max(dsid)))';         % different A2 for each dataset
    tau2 = params((3+2*max(dsid)):(2+3*max(dsid)))';   % different tau2 for each dataset
    
    F = A0(dsid)+A1*exp(-T/tau1)+A2(dsid).*exp(-T./tau2(dsid));
   
    