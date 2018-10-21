
function [c, ceq, dc, dceq] = constFun(X, aircraft, N)
   
    n_coeffs = 2*N+1;
    coeffs_x = X(1:n_coeffs,1);
    coeffs_y = X(n_coeffs+1:2*n_coeffs,1);
    coeffs_z = X(2*n_coeffs+1:3*n_coeffs,1);
    coeffs = [coeffs_x,coeffs_y, coeffs_z];
    VR = X(end-1,1); tf = X(end,1);
    
    M = 7*N; c = zeros(13*(M+1),1);
    t = linspace(0, tf, M);
    for i = 1:M
        j = (i-1)*13 + 1;
        % get flat outputs at time t
        sigma = aircraft.get_traj(t(i), tf, coeffs, N); 
        z = sigma(3); 
        % get states at time t`
        aircraft = aircraft.get_xu(sigma, VR);
        % constrain states and controls
        V = aircraft.x(1); Cl = aircraft.u(1); nu = aircraft.u(2);
        CT = aircraft.u(3);
        % Cl constraints
        c(j,1) = aircraft.limits(1) - Cl;
        c(j+1,1) = Cl - aircraft.limits(2);
        % V contraints
        c(j+2,1) = aircraft.limits(3) - V;
        c(j+3,1) = V - aircraft.limits(4);
        % nu constraints
        c(j+4,1) = aircraft.limits(5) - nu;
        c(j+5,1) = nu - aircraft.limits(6);
        % CT constraints
        c(j+6,1) = aircraft.limits(7) - CT;
        c(j+7,1) = CT - aircraft.limits(8);
        % hmin constraints
        c(j+8,1)  = z + 0.5*aircraft.b*sin(nu);
        c(j+9,1) = z - 0.5*aircraft.b*sin(nu);
        c(j+10,1) = -z - 100;
        % gamma constraints
        c(j+11,1) = aircraft.limits(9) - aircraft.x(2);
        c(j+12,1) = aircraft.x(2) - aircraft.limits(10);
        
    end
    
    % equality constraints to ensure periodicity
    sigma_0 = aircraft.get_traj(0, tf, coeffs, N);
    sigma_f = aircraft.get_traj(tf, tf, coeffs, N);
    ceq(1,1) = sigma_0(1) - sigma_f(1);
    ceq(2,1) = sigma_0(2) - sigma_f(2);
    ceq(3,1) = sigma_0(3) - sigma_f(3);
    
    if nargout > 2 % gradient of the constraints
      dc = [];
      dceq = [];
    end        
end

