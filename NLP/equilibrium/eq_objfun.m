function [f, g] = eq_objfun(x, p)
    fvec = sys_model(0, x, p);
    f = fvec'*fvec;
    
    if nargout > 1
        J = get_jac(x, p);
        g = 2*J'*fvec;
    end
end