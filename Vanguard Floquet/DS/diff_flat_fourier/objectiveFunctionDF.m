function f = objectiveFunctionDF(X)
    N = (length(X)-2)/3;
    f = X(3*N+2);
end