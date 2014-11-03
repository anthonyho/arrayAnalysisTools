function F = DEF10302014(params,xdata)
% For fitting 10/30/2014, scf11

    %%
    
    i = 6;
    
    allFixedParams = [ ...
    0.1599    0.2288    0.1446    0.1297    0.1255    0.1461
    0.2829    0.3833    0.3165    0.4219    0.4115    0.3066
    0.5607    0.4140    0.5454    0.4558    0.4687    0.5504];

    fixedParams = allFixedParams(:,i);

    %%
    
    F = exp(-xdata/params(1)) .* ( ...
        fixedParams(3)*exp(-xdata/params(2)) + fixedParams(2) ) ...
        + fixedParams(1);

end