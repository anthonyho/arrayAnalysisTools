function F = doubleExponential4Distinct(x,xdata)
    F = (1-x(1)-x(3))+x(1)*exp(-xdata/(x(2)^2+x(4)))+x(3)*exp(-xdata/x(4));