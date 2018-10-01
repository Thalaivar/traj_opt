function f = objFun(x, p, A)
    N = p(1);
    eig_val = x(1,1);
    eig_vec = gen_vec(x(2:end,1), N);   
    v1 = (A - eig_val*eye(N))*eig_vec;
    v2 = (norm(eig_vec) - 1)^2;
    f = v1'*v1 + v2^2;
end