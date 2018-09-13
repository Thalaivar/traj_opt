function [f, g] = eq_objfun(x)
    fvec = sys_model(0, x);
    f = fvec'*fvec;
    
    if nargout > 1
        J = get_jac(x);
        g = 2*J'*fvec;
    end
end