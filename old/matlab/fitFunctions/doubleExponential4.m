function F = doubleExponential4(x,xdata)
    F = (1-x(1)-x(3))+x(1)*exp(-xdata/x(2))+x(3)*exp(-xdata/x(4));