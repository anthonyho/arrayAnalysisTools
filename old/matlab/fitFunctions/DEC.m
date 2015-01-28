function F = DEC(x,xdata)
    F = x(1)+x(2)*exp(-xdata/x(3))+(1-x(1)-x(2))*exp(-xdata/x(4));
end