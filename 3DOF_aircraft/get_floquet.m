function [eigvec, eigval, sol] = get_floquet(aircraft, N, M, xguess)
    
    % dimension of eigenvec
    n = 3;  
    % dec vec is [re(u1), re(u2), .... , re(uM), im(u1), ... , re(lam), im(lam)]
    if isempty(xguess)
        xguess = [(1/sqrt(6))*ones(2*M*n,1);0.0;0.0];
    end
    
    lb = ones(2*M*n+2,1);
    ub = ones(2*M*n+2,1);
    
    % bounds for eig_vec components
    lb(1:2*M*n,1) = -100*lb(1:2*M*n,1); 
    ub(1:2*M*n,1) = 100*ub(1:2*M*n,1);
    % bounds for eigenvalue
    lb(end-1,1) = -Inf; ub(end-1,1) = Inf;
    lb(end,1) = -Inf; ub(end,1) = Inf;
    
    options = optimoptions('fmincon', 'Display', 'Iter', 'Algorithm', 'sqp', 'MaxFunctionEvaluations', 1000000, 'ConstraintTolerance', 1e-5);
    sol = fmincon(@(x) objfun(x, N,2), xguess, [], [], [], [], lb, ub, @(x) constFun(x, aircraft, N, M, true), options);
    
    % retrieve eigenvectors at each point
    eigvec_comp = zeros(M,6);
    for i = 1:6
        j = (i-1)*M + 1;
        eigvec_comp(:,i) = sol(j:j+M-1,1)';
    end
    alpha = zeros(3,M); beta = zeros(3,M);
    for i = 1:M
        alpha(:,i) = eigvec_comp(i,1:3)';
        beta(:,i) = eigvec_comp(i,4:6)';
    end
    eigvec = zeros(3,M);
    for i = 1:M
        eigvec(:,i) = complex(alpha(:,i), beta(:,i));
    end
    
    % retrieve eigenvalue
    eigval = complex(sol(end-1), sol(end));
    
    
    
end