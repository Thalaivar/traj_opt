function [a, eta, tf, VR] = generate_traj(limits, model_par, N)
    % limits is of the form:
    %       limits = [Clmax, Vmax, nu_min, nu_max, Tmin, Tmax, hmin]    
    
    % model_par is of the form:
    %       model_par = [m, rho, S, g, b, Cd0, Cd1, Cd2] 
    
    % dec vec = [ah_0, ..., ah_N, ax_0, ...,ax_N, ay_0, ..., ay_N, .., VR, tf]
    n_coeffs = N+1; n_phase_angles = N;
    x0 = zeros(3*(n_coeffs+n_phase_angles) + 2,1);
    % setting initial guess params for coeffs
    x0(1,1) = 20; x0(2,1) = -20; 
    x0(n_coeffs+1,1) = -20; x0(n_coeffs+2,1) = 20; 
    x0(2*n_coeffs+1,1) = 0; x0(2*n_coeffs+2,1) = -20*(2^0.5);
    % setting initial guess params for phase angles
    x0(3*n_coeffs+1,1) = pi/2;
    x0(3*n_coeffs+n_phase_angles+1,1) = pi/2;
    x0(3*n_coeffs+2*n_phase_angles+1,1) = 0;
    % setting initial guess params for VR, tf
    x0(end-1,1) = 12; x0(end,1) = 10;
    lb = ones(3*(n_coeffs+n_phase_angles) + 2,1,1); ub = ones(3*(n_coeffs+n_phase_angles) + 2,1,1);
    % coeffs of trajectory
    lb(1:n_coeffs*3,1) = -500*lb(1:n_coeffs*3,1); ub(1:n_coeffs*3,1) = 500*ub(1:n_coeffs*3,1);
    % phase angles
    lb(n_coeffs*3+1:3*(n_coeffs+n_phase_angles),1) = 0*lb(n_coeffs*3+1:3*(n_coeffs+n_phase_angles),1); 
    ub(n_coeffs*3+1:3*(n_coeffs+n_phase_angles),1) = 2*pi*ub(n_coeffs*3+1:3*(n_coeffs+n_phase_angles),1);
    % VR
    lb(end-1,1) = 0; ub(end-1,1) = 150;
    % tf
    lb(end,1) = 0; ub(end,1) = 200;
    
    [c, ceq] = constFun(x0, limits, model_par, N);
    
    options = optimoptions('fmincon', 'Display', 'Iter', 'Algorithm', 'sqp', 'MaxFunctionEvaluations', 10000, 'StepTolerance', 1e-10, 'ConstraintTolerance', 1);
    [x, fval] = fmincon(@(x) objfun(x), x0, [], [], [], [], lb, ub, @(x) constFun(x, limits, model_par, N), options);
    
    n1 = n_coeffs; n2 = n_phase_angles;
    a = [x(1:n1,1),x(n1+1:2*n1,1),x(2*n1+1:3*n1,1)];
    eta = [x(3*n1+1:3*n1+n2,1), x(3*n1+n2+1:3*n1+2*n2,1), x(3*n1+2*n2+1:3*(n1+n2),1)];    
    VR = x(end-1,1); tf = x(end,1);
end