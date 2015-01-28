function F = singleExponential(x,xdata)
    F = x(1)+x(2)*exp(-xdata/x(3));