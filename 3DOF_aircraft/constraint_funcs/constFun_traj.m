% to be used when optimising for periodic trajectories ONLY
% NOTE: 
%   aircraft needs N and params set (and VR for stability optimisation)
function [c , ceq] = constFun_traj(X, aircraft, type, zmax_nom)
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
       
       M = 1000; % no. of points at which dynamic constraints are imposed

       c = zeros(16*M,1);
       t = linspace(0, tf, M);
       for i = 1:M
           j = (i-1)*16 + 1;
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
           c(j+2,1) = 10 - V;      % 12.9372
           c(j+3,1) = V - 80; % 31.1240
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
           c(j+10,1) = -pi/4 - gamma; % -0.5549
           c(j+11,1) = gamma - pi/4; % 0.4011
           % x and y constraints
           c(j+12,1) = -500 - x; % -500
           c(j+13,1) = x - 500;  % 500
           c(j+14,1) = -500 - y; % -500 
           c(j+15,1) = y - 500;  % 500
       end
        
        if(~isempty(zmax_nom))
         % matching the zmin for nominal and unstable solution
           t = linspace(0, tf, 1000);
           zmax = -Inf;
           for i = 1:length(t)
               sig = get_traj(t(i), tf, coeffs, N);
               if sig(3) >  zmax, zmax = sig(3); end
           end
           ceq = zmax - zmax_nom;
        else, ceq = [];
        end 
       
end