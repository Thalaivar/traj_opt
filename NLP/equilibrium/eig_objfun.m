function f = eig_objfun(x, J, ch)
    N = size(J); N = N(1);
    eig_val = x(1:2*N,1);
    eigvec_comp = x(2*N+1:end,1); % components of eigenvectors
    eig_vec = gen_vec(eigvec_comp, N); % eig_vals = [v1, v2, ..., vn]
    f = 0;
    for i = 1:2:N
        %v1 = (J - complex(eig_val(i,1), eig_val(i+1,1))*eye(N))*eig_vec(:,i);
        v2 = det(J - complex(eig_val(i), eig_val(i+1,1))*eye(N));
        f = f + v2*conj(v2);
    end
end