function [a, eta, tf, VR] = generate_traj(limits, model_par, N)
    % limits is of the form:
    %       limits = [Clmax, Vmax, nu_min, nu_max, Tmin, Tmax, hmin]    
    
    % model_par is of the form:
    %       model_par = [m, rho, S, g, b, Cd0, Cd1, Cd2] 
    
    % dec vec = [ah_0, ..., ah_N, ax_0, ...,ax_N, ay_0, ..., ay_N, .., VR, tf]
    
    x0 = ones((N+1)*6 + 2,1); x0(end,1) = 1; x0(end-1,1) = 1;
    lb = ones((N+1)*6 + 2,1); ub = ones((N+1)*6 + 2,1);
    % coeffs of trajectory
    lb(1:(N+1)*3,1) = -500*lb(1:(N+1)*3,1); ub(1:(N+1)*3,1) = 500*ub(1:(N+1)*3,1);
    % phase angles
    lb((N+1)*3+1:(N+1)*6,1) = 0*lb((N+1)*3+1:(N+1)*6,1); ub((N+1)*3+1:(N+1)*6,1) = 2*pi*ub((N+1)*3+1:(N+1)*6,1);
    % VR
    lb(end-1,1) = 0; ub(end-1,1) = 150;
    % tf
    lb(end,1) = 0; ub(end,1) = 200;
    
    [c, ceq] = constFun(x0, limits, model_par, N)
    
    options = optimoptions('fmincon', 'Display', 'Iter', 'Algorithm', 'sqp', 'MaxFunctionEvaluations', 10000);
    [x, fval] = fmincon(@(x) objfun(x), x0, [], [], [], [], lb, ub, @(x) constFun(x, limits, model_par, N), options)
    
    a = [x(1:N+1,1),x(N+2:2*(N+1),1),x(2*N+3,3*(N+1),1)];
    eta = [x(3*N+4:4*(N+1),1), x(4*N+5:5*(N+1),1), x(5*N+6:6*(N+1),1)];    
    VR = x(end-1,1); tf = x(end,1);
end