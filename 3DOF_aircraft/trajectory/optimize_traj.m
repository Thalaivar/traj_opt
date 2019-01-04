% to generate optimal periodic trajectories
% aircraft should have N, params set
function [aircraft, sol] = optimize_traj(aircraft, x0, p)    
    % dec vec = [ax0, ... , ax1_N,ax2_1, ... , ax2_N, ay0, ... , ay2_N, az0, ... ,az2_N, VR, tf]
    % a$1_i corresponds to teh coefficient of the i'th cosine harmonic
    % a$2_i corresponds to teh coefficient of the i'th sine harmonic
    xguess = x0;
    aircraft.p = p; N = aircraft.N;
    
    % if x0 was empty, it means we are running for first time
    if isempty(xguess)
        xguess = get_init_guess('circle', N);
        "INITIAL GUESS WAS NOT PROVIDED!"
    end
    
    lb = ones(3*(2*N+1)+2,1); ub = ones(3*(2*N+1)+2,1);
    % bounds on coefficients
    lb(1:3*(2*N+1),1) = -500*lb(1:3*(2*N+1),1);
    ub(1:3*(2*N+1),1) = 500*ub(1:3*(2*N+1),1);
    % bounds on VR and tf
    lb(end-1,1) = 0; ub(end,1) = 150;
    lb(end,1) = 0; ub(end,1) = 200;
    
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point', 'MaxFunctionEvaluations', 2000000, 'MaxIterations', 10000);
    sol = fmincon(@(x) objfun(x, aircraft, 'traj'), xguess, [], [], [], [], lb, ub, @(x) constFun_traj(x, aircraft, 'traj'), options);
    
    n = (2*N+1);
    aircraft.coeffs = [sol(1:n,1), sol(n+1:2*n,1), sol(2*n+1:3*n,1)];
    aircraft.VR = sol(3*n+1,1); aircraft.tf = sol(3*n+2,1);
end