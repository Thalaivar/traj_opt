function [a, eta, tf, VR, eig_vec, eig_val, x] = get_floquet(limits, model_par, xguess, N, M)
    % model_par is of the form:
    %       model_par = [m, rho, S, g,Cd0, Cd1, Cd2, b]
    
     % limits is of the form:
    %       limits = [Clmax, Vmax, nu_min, nu_max, Tmin, Tmax, hmin]    
    
    
    % n1 : n_coeffs; n2 : n_phase_angles; n3 : n_eigvec_comp;
    n1 = N+1; n2 = N;
    n3 = 6*M;
    
    
    x0 = [xguess;ones(n3,1);0];
    lb = ones(3*(n1+n2)+2+n3+1,1);
    ub = ones(3*(n1+n2)+2+n3+1,1);
    % bounds for coeffs
    lb(1:3*n1,1) = -500*lb(1:3*n1,1); ub(1:3*n1,1) = 500*ub(1:3*n1,1);
    % bounds for phase angles
    lb(3*n1+1:3*(n1+n2),1) = 0*lb(3*n1+1:3*(n1+n2),1); 
    ub(3*n1+1:3*(n1+n2),1) = 2*pi*ub(3*n1+1:3*(n1+n2),1); 
    % bounds for VR, tf  ]                                                
    lb(3*(n1+n2)+1,1) = 0*lb(3*(n1+n2)+1,1); ub(3*(n1+n2)+1,1) = 100*ub(3*(n1+n2)+1,1);
    lb(3*(n1+n2)+2,1) = 0*lb(3*(n1+n2)+2,1); ub(3*(n1+n2)+2,1) = 200*ub(3*(n1+n2)+2,1);
    % bounds for eig_vec components
    lb(3*(n1+n2)+2+1:3*(n1+n2)+2+n3,1) = -100*lb(3*(n1+n2)+2+1:3*(n1+n2)+2+n3,1); 
    ub(3*(n1+n2)+2+1:3*(n1+n2)+2+n3,1) = 100*ub(3*(n1+n2)+2+1:3*(n1+n2)+2+n3,1);
    % bounds for eigenvalue
    lb(end,1) = -10; ub(end,1) = 10;
    
    options = optimoptions('fmincon', 'Display', 'Iter', 'Algorithm', 'sqp', 'MaxFunctionEvaluations', 1000000, 'StepTolerance', 1e-10, 'MaxIterations', 10000);
    x = fmincon(@(x) objfun2(x,N), x0, [], [], [], [], lb, ub, @(x) constFun_floq(x, limits, model_par, N, M), options);
    
    a = [x(1:n1,1),x(n1+1:2*n1,1),x(2*n1+1:3*n1,1)];
    eta = [x(3*n1+1:3*n1+n2,1), x(3*n1+n2+1:3*n1+2*n2,1), x(3*n1+2*n2+1:3*(n1+n2),1)];    
    VR = x(3*(n1+n2)+1,1); tf = x(3*(n1+n2)+2,1);
    eig_vec = x(3*(n1+n2)+2+1:3*(n1+n2)+2+n3,1);
    eig_val = x(end,1);
end