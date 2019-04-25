function f = objfun_colloc(X)
    N = (length(X)-2)/8;
    
    f = -X(8*N+1,1);
end