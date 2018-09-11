function f = objfun(x, N)
    f1 = x(end,1)^2;
    f2 = (x(1,1) - x(N,1))^2 + (x(N+1,1) - x(2*N,1) + 2*pi)^2;
    %f = f1;
    f = f2;
    %f = f1+f2;
end