function [eig_vec, eig_val, sol] = get_floquet(aircraft, N, M, xguess)
    
    % dimension of eigenvec
    n = 6;  
    % dec vec is [re(u1), re(u2), .... , re(uM), im(u1), ... , re(lam), im(lam)]
    if isempty(xguess)
        xguess = [(1/sqrt(12))*ones(2*M*n,1);0.7;0.9852];
    end
    
    lb = ones(2*M*n+2,1);
    ub = ones(2*M*n+2,1);
    
    % bounds for eig_vec components
    lb(1:2*M*n,1) = -1*lb(1:2*M*n,1); 
    ub(1:2*M*n,1) = 1*ub(1:2*M*n,1);
    % bounds for eigenvalue
    lb(end-1,1) = -10; ub(end-1,1) = 10;
    lb(end,1) = -10; ub(end-1,1) = 10;
    
    options = optimoptions('fmincon', 'Display', 'Iter', 'Algorithm', 'interior-point', 'MaxFunctionEvaluations', 100000);
    sol = fmincon(@(x) objfun(x, N,2), xguess, [], [], [], [], lb, ub, @(x) constFun(x, aircraft, N, M, true), options);
    
    eigvec_comp_real = sol(1:M*n,1); eigvec_comp_img = sol(M*n+1:2*M*n,1);
    eig_val = [sol(end-1,1), sol(end,1)];
    eig_vec = zeros(6,M); 
       for i = 1:M
           j = (i-1)*6 + 1;
           eig_vec(:,i) = complex(eigvec_comp_real(j:j+5,1), eigvec_comp_img(j:j+5,1));
       end
end