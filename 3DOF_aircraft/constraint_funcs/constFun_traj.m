% to be used when optimising for periodic trajectories ONLY
% NOTE: 
%   aircraft needs N and params set
function [c , ceq] = constFun_traj(X, aircraft)
       N = aircraft.N;
       n_coeffs = 2*N+1;
       coeffs_x = X(1:n_coeffs,1);
       coeffs_y = X(n_coeffs+1:2*n_coeffs,1);
       coeffs_z = X(2*n_coeffs+1:3*n_coeffs,1);
       coeffs = [coeffs_x,coeffs_y, coeffs_z];
       VR = X(3*n_coeffs+1,1); tf = X(3*n_coeffs+2,1);
       aircraft.tf = tf; aircraft.VR = VR;
       M = 75; % no. of points at which dynamic constraints are imposed

       c = zeros(18*M,1);
       t = linspace(0, tf, M);
       for i = 1:M
           j = (i-1)*18 + 1;
           % get flat outputs at time t
           sigma = get_traj(t(i), tf, coeffs, N);
           x = sigma(1); y = sigma(2); z = sigma (3);
           % get states at time t`
           aircraft = aircraft.get_xu(sigma);
           % constrain states and controls
           V = aircraft.x(1); Cl = aircraft.u(1); nu = aircraft.u(2);
           CT = aircraft.u(3); gamma = aircraft.x(2);
           % Cl constraints
           c(j,1) = -0.2 - Cl;
           c(j+1,1) = Cl - 1.17;
           % V contraints
           c(j+2,1) = 10 - V;
           c(j+3,1) = V - 80;
           % nu constraints
           c(j+4,1) = -pi/3 - nu;
           c(j+5,1) = nu - pi/3;
           % CT constraints
           c(j+6,1) = -1e-4 - CT;
           c(j+7,1) = CT - 1e-4;
           % hmin constraints
           c(j+8,1)  = z + 0.5*aircraft.b*sin(nu);
           c(j+9,1) = z - 0.5*aircraft.b*sin(nu);
           % gamma constraints
           c(j+10,1) = -pi/4 - gamma;
           c(j+11,1) = gamma - pi/4;
           % x and y constraints
           c(j+12,1) = -500 - x;
           c(j+13,1) = x - 500;
           c(j+14,1) = -500 - y;
           c(j+15,1) = y - 500;
           % z constraints
           c(j+16,1) = z + 0.5*aircraft.b*abs(sin(nu));
           c(j+17,1) = -100 - z;
            
       end

       ceq = [];
end