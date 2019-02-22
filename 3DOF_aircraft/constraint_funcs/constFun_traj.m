% to be used when optimising for periodic trajectories ONLY
% NOTE: 
%   aircraft needs N and params set (and VR for stability optimisation)
function [c , ceq] = constFun_traj(X, aircraft, type, ref_ac)
       N = aircraft.N;
       n_coeffs = 2*N+1;
       coeffs_x = X(1:n_coeffs,1);
       coeffs_y = X(n_coeffs+1:2*n_coeffs,1);
       coeffs_z = X(2*n_coeffs+1:3*n_coeffs,1);
       coeffs = [coeffs_x,coeffs_y, coeffs_z];
        
       if strcmp(type, 'traj')
           VR = X(3*n_coeffs+1,1); tf = X(3*n_coeffs+2,1);
           aircraft.tf = tf; aircraft.VR = VR;
       elseif strcmp(type, 'stability')
           tf = X(3*n_coeffs+1,1); aircraft.tf = tf;
       end
       
       M = 50; % no. of points at which dynamic constraints are imposed

       c = zeros(18*M,1);
       t = linspace(0, tf, M);
       t_nom = linspace(0,
       for i = 1:M
           j = (i-1)*18 + 1;
           % get flat outputs at time t
           sigma = get_traj(t(i), tf, coeffs, N);
           sig_nom = get_traj(t(i), ref_ac.tf, ref_ac.coeffs, ref_ac.N);
           
           x = sigma(1); y = sigma(2); z = sigma (3);
           x_nom = sig_nom(1); y_nom = sig_nom(2); z_nom = sig_nom(3);
           
           % get states at time t`
           aircraft = aircraft.get_xu(sigma);
           % constrain states and controls
           V = aircraft.x(1); Cl = aircraft.u(1); nu = aircraft.u(2);
           CT = aircraft.u(3); gamma = aircraft.x(2);
           
           ref_ac = ref_ac.get_xu(sig_nom);
           V_nom = ref_ac.x(1); gamma_nom = ref_ac.x(3);
           
           % Cl constraints
           c(j,1) = -0.2 - Cl;
           c(j+1,1) = Cl - 1.17;
           % V contraints
           c(j+2,1) = 0.4*V_nom - V;      % 10
           c(j+3,1) = V - 1.6*V_nom; % 80
           % nu constraints
           c(j+4,1) = -pi/3 - nu;
           c(j+5,1) = nu - pi/3;
           % CT constraints
           c(j+6,1) = -1e-4 - CT;
           c(j+7,1) = CT - 1e-4;
           % hmin constraints
           hmin = 0.15;
           c(j+8,1)  = z + 0.5*aircraft.b*sin(nu) + hmin;
           c(j+9,1) = z - 0.5*aircraft.b*sin(nu) + hmin;
           % gamma constraints
           c(j+10,1) = 0.4*gamma_nom - gamma; % -pi/4
           c(j+11,1) = gamma - 1.6*gamma_nom; % pi/4
           % x and y constraints
           c(j+12,1) = 0.4*x_nom- x; % -500
           c(j+13,1) = x - 1.6*x_nom;  % 500
           c(j+14,1) = 0.4*y_nom - y; % -500 
           c(j+15,1) = y - 1.6*y_nom;  % 500
           % z constraints
           c(j+16,1) = 0.4*z_nom - z;
           c(j+17,1) = z - 1.6*z_nom;
       end

       ceq = [];
end