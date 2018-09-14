clear all
clc

p = [2/3, 4/3, 1, 1];

% get equilibrium point
xguess = [0;0];
xeq = find_equilibrium(xguess, p);

J = get_jac(xeq, p)
N = size(J); N = N(1);
% get eigenvalues
% state vector is x = [ λ1, λ2, ... , λn, v1, v2, ... , vn]'
n_eigvals = N; n_eigvecs = N;
x0 = ones(n_eigvals + n_eigvecs*N,1);
for i = 1:n_eigvals
    if mod(i, 2) == 0
        x0(i) = -1;
    end
end
x0(n_eigvals+1:end,1) = 2^(-0.5)*x0(n_eigvals+1:end,1);
options = optimoptions('fmincon', 'Display', 'iter', 'MaxFunctionEvaluations', 200000, 'MaxIterations', 10000);
[x, fval] = fmincon(@(x) eig_objfun(x, J), x0, [], [], [], [], [], [], @(x) constFun(x, J));

eig_val = x(1:n_eigvals,1);
eigvec_comp = x(n_eigvals+1:end,1);
eig_vec = gen_vec(eigvec_comp, N);
    

