function f = objfun(x, N, type)
    if(type == 1)
        f = x(3*(2*N+1)+1,1);
    end
    
    if(type == 2)
        f = -x(end-1,1);
    end
        
end