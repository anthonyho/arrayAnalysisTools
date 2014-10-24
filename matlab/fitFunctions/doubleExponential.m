function F = doubleExponential(x,xdata)
    F = x(1)+x(2)*exp(-xdata/x(3))+x(4)*exp(-xdata/x(5));