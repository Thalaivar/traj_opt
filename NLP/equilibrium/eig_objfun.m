function f = eig_objfun(x, J)
    N = size(J);
    N = N(1);
    eig_val = x(1:N,1);
    eigvec_comp = x(N+1:end,1);
    fval = 0;
    j = 1;
    for i = 1:N
        j = 2*i - 1;
        eig_vec = [eigvec_comp(j); eigvec_comp(j+1)];
        v = (J - eig_val(i,1)*eye(N))*eig_vec;
        fval = fval + v'*v;
    end
    f = fval'*fval;
    
end