function F = globalDoubleExponential4(params,X)
    T = X(:,1);        % time
    dsid = X(:,2);     % dataset id
    A1 = params(1:max(dsid))';          % different A0 for each datasets
    %A1 = params(1+max(dsid));           % same A1 for all datasets
    tau1 = params(1+max(dsid));         % same tau1 for all datasets
    A2 = params((2+max(dsid)):(1+2*max(dsid)))';         % different A2 for each dataset
    tau2 = params((2+2*max(dsid)):(1+3*max(dsid)))';   % different tau2 for each dataset
    
    F = (1-A1(dsid)-A2(dsid))+A1(dsid).*exp(-T/tau1)+A2(dsid).*exp(-T./tau2(dsid));
   
    