%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to calculate states w/o flatness                      %
% INPUT :                                                        %
%       t        - time                                          %
%       X        - state vector: [V, chi, gamma, x, y, z]'       %
%       solution - struct with fields 'VR', 'tf', 'coeffs', 'N'  %
% OUTPUT:                                                        %
%       ydot - derivative of states                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ydot = non_flat_model(t, X, solution)
    % model constants
    m = 4.5; rho = 1.225; S = 0.473;
    g = 9.806; Cd0 = 0.0173; Cd1 = -0.0337;
    Cd2 = 0.0517;
    
    VR = solution.VR; tf = solution.tf; 
    coeffs = solution.coeffs; N = solution.N;
    
    V = X(1,1); gamma = X(3,1); chi = X(2,1);
    z = X(6,1);  
    
    % get control inputs from time history
    sigma = get_traj(t, tf, coeffs, N);
    [~, u] = get_xu(sigma, VR);
    Cl = u(1); nu = u(2); CT = u(3);

    Cd = Cd0 + Cd1*Cl + Cd2*Cl^2;
    D = 0.5*rho*S*V^2*Cd;
    L = 0.5*rho*S*V^2*Cl;
    T = 0.5*rho*S*V^2*CT;

     % wind model
    p_exp = 1;
    Wx = VR*(-z)^p_exp;
    Wxz = (p_exp*VR)*((-z)^p_exp)/z;

    zdot = -V*sin(gamma); xdot = V*cos(chi)*cos(gamma) + Wx; ydot = V*sin(chi)*cos(gamma);
    Vdot = (-D/m) - g*sin(gamma) - Wxz*zdot*cos(chi)*cos(gamma) + (T/m);
    chidot = (L*sin(nu)/(m*V*cos(gamma))) + Wxz*zdot*sin(chi)/(V*cos(gamma));
    gammadot = (L*cos(nu)/(m*V)) - g*cos(gamma)/V + Wxz*zdot*cos(chi)*sin(gamma)/V;

    ydot = [Vdot; chidot; gammadot; xdot; ydot; zdot];
 end
