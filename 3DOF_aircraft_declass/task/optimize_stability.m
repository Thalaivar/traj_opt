% to generate periodic traj with optimal stability
% aircraft should have VR, N and params set
function [sol, X] = optimize_stability(x0, N)    
    % dec vec = [ax0, ... , ax1_N,ax2_1, ... , ax2_N, ay0, ... , ay2_N, az0, ... ,az2_N, tf]
    % a$1_i corresponds to teh coefficient of the i'th cosine harmonic
    % a$2_i corresponds to teh coefficient of the i'th sine harmonic

    % if x0 was empty, it means we are running for first time
    if isempty(x0)
        xguess = get_init_guess('circle', N);
    else
        xguess = x0;
    end
    solution.VR = xguess(end-1,1); solution.N = N; 
    solution.tf = xguess(end,1);
    n = (2*N+1);
    solution.coeffs = [xguess(1:n,1), xguess(n+1:2*n,1), xguess(2*n+1:3*n,1)];  
    
    lb = ones(3*(2*N+1)+2,1); ub = ones(3*(2*N+1)+2,1);
    % bounds on coefficients
    lb(1:3*(2*N+1),1) = -500*lb(1:3*(2*N+1),1);
    ub(1:3*(2*N+1),1) = 500*ub(1:3*(2*N+1),1);
    % bounds on VR and tf
    lb(end-1,1) = 0; ub(end-1,1) = 100;
    lb(end,1) = 0; ub(end,1) = 150;
    
    options = optimoptions('fmincon', 'Display', 'Iter', 'Algorithm', 'interior-point1', 'UseParallel', true);
    options.MaxFunctionEvaluations = 100000;
    options.StepTolerance = 1e-15;
    options.MaxIterations = 10000;
    
    solution.objfun_type = 'floq_new';
    solution.constFun_type = 'traj';
    
    X = fmincon(@(x) objfun(x, solution), xguess, [], [], [], [], lb, ub, @(x) constFun_traj(x, solution), options);
    
    n = (2*N+1);
    sol.coeffs = [X(1:n,1), X(n+1:2*n,1), X(2*n+1:3*n,1)];
    sol.tf = X(3*n+2,1); sol.VR = X(3*n+1,1);
end