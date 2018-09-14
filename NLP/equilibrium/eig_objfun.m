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
        v1 = (J - eig_val(i,1)*eye(N))*eig_vec;
        v2 = det(J - eig_val(i)*eye(N));
        fval = fval + v1'*v1 + v2^2;
    end
    f = fval'*fval;
    
end