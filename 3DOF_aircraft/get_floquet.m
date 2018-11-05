function [eig_vec, eig_val, sol] = get_floquet(aircraft, N, M, xguess)
    
    % dimension of eigenvec
    n = 3;  
    % dec vec is [re(u1), re(u2), .... , re(uM), im(u1), ... , re(lam), im(lam)]
    if isempty(xguess)
        xguess = [(1/sqrt(6))*ones(2*M*n,1);0.01;0.0];
    end
    
    lb = ones(2*M*n+2,1);
    ub = ones(2*M*n+2,1);
    
    % bounds for eig_vec components
    lb(1:2*M*n,1) = -100*lb(1:2*M*n,1); 
    ub(1:2*M*n,1) = 100*ub(1:2*M*n,1);
    % bounds for eigenvalue
    lb(end-1,1) = -100; ub(end-1,1) = Inf;
    lb(end,1) = -100; ub(end,1) = Inf;
    
    options = optimoptions('fmincon', 'Display', 'Iter', 'Algorithm', 'sqp', 'MaxFunctionEvaluations', 1000000);
    sol = fmincon(@(x) objfun(x, N,2), xguess, [], [], [], [], lb, ub, @(x) constFun(x, aircraft, N, M, true), options);
    
    eigvec_comp_real = sol(1:M*n,1); eigvec_comp_img = sol(M*n+1:2*M*n,1);
    eig_val = [sol(end-1,1), sol(end,1)];
    eig_vec = zeros(3,M); 
       for i = 1:M
           j = (i-1)*3 + 1;
           eig_vec(:,i) = complex(eigvec_comp_real(j:j+2,1), eigvec_comp_img(j:j+2,1));
       end
end