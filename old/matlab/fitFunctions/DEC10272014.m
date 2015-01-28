function F = DEC10272014(params,xdata)
% For fitting 10/27/2014, scf5

    %%

    F = exp(-xdata/params(3)) .* ( ...
        (params(2)+params(4))*exp(-xdata/params(5)) + params(2) ) ...
        + params(1) - 2*params(2) - params(4);

end