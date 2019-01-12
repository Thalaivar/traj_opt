%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to calculate states and control inputs from flatness  %
% INPUT :                                                        %
%       sigma - flat outputs                                     %
%       VR    - reference wind speed                             %
% OUTPUT:                                                        %
%       x - wind axes states : [V, gamma, chi]                   %
%       u - control inputs   : [Cl, nu, CT]                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, u] = get_xu(sigma, VR)
    m = 4.5; rho = 1.225; S = 0.473;
    g = 9.806; Cd0 = 0.0173; Cd1 = -0.0337;
    Cd2 = 0.0517; b = 3;
    z = sigma(3);
    zdot = sigma(6); xdot = sigma(4); ydot = sigma(5); 
    zddot = sigma(9); xddot = sigma(7); yddot = sigma(8);
    
    % wind model
    p_exp = 1;
    Wx = VR*(-z)^p_exp;
    Wxz = (p_exp*VR)*((-z)^p_exp)/z;

    % non flat outputs
    V = ((xdot - Wx)^2 + ydot^2 + zdot^2)^0.5;
    Vdot = (xdot*xddot - xdot*zdot*Wxz - xddot*Wx + Wx*Wxz*zdot + ydot*yddot + zdot*zddot)/V;
    gamma = -asin(zdot/V);
    gammadot = (zdot*Vdot - V*zddot)/(V*(V^2 - zdot^2)^0.5);
    chi = atan2(ydot,(xdot - Wx));
    chidot = (xdot*yddot - yddot*Wx - ydot*xddot + ydot*zdot*Wxz)/(ydot^2 + xdot^2 + Wx^2 - 2*xdot*Wx);
    nu = atan((V*cos(gamma)*chidot - Wxz*zdot*sin(chi))/(V*gammadot + g*cos(gamma) - Wxz*cos(chi)*sin(gamma)*zdot));
    Cl = (m*V*cos(gamma)*chidot - m*Wxz*zdot*sin(chi))/(0.5*rho*S*sin(nu)*V^2);
%             if(nu == 0 && (chidot == 0 || zdot == 0 || chi == 0))
%                 Cl = 0;
%             end

    % aerodynamic forces
    Cd = Cd0 + Cd1*Cl + Cd2*Cl^2;
    D = 0.5*rho*S*V^2*Cd;
    T = m*Vdot + D + m*g*sin(gamma) + m*Wxz*zdot*cos(gamma)*cos(chi);
    CT = T/(0.5*rho*S*V^2);
    x = [V, gamma, chi]; u = [Cl, nu, CT];
end