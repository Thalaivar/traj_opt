function f = eig_objfun(x, J, ch)
    N = size(J); N = N(1);
    eig_val = x(1:N,1);
    eigvec_comp = x(N+1:end,1); % components of eigenvectors
    eig_vec = gen_vec(eigvec_comp, N); % eig_vals = [v1, v2, ..., vn]
    f = 0;
    for i = 1:length(eig_val)
        v1 = (J - eig_val(i,1)*eye(N))*eig_vec(:,i);
        v2 = det(J - eig_val(i)*eye(N));
        f = f + v1'*v1 + v2^2;
    end
end