% to generate optimal periodic trajectories
% aircraft should have N, params set
function [sol, X] = optimize_traj(x0, N)    
    % dec vec = [ax0, ... , ax1_N,ax2_1, ... , ax2_N, ay0, ... , ay2_N, az0, ... ,az2_N, VR, tf]
    % a$1_i corresponds to teh coefficient of the i'th cosine harmonic
    % a$2_i corresponds to teh coefficient of the i'th sine harmonic    
    
    % if x0 was empty, it means we are running for first time
    if isempty(x0)
        xguess = get_init_guess('circle', N);
        "INITIAL GUESS WAS NOT PROVIDED!"
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
    lb(end-1,1) = 0; ub(end,1) = 150;
    lb(end,1) = 0; ub(end,1) = 200;
    
    % type of optimisation
    solution.constFun_type = 'traj';
    solution.objfun_type   = 'traj';
    
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point', 'MaxFunctionEvaluations', 2000000, 'MaxIterations', 10000);
    X = fmincon(@(x) objfun(x, solution), xguess, [], [], [], [], lb, ub, @(x) constFun_traj(x, solution), options);
  
    sol.coeffs = [X(1:n,1), X(n+1:2*n,1), X(2*n+1:3*n,1)];
    sol.VR = X(3*n+1,1); sol.tf = X(3*n+2,1);
end