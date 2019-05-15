function [state, control] = get_traj(t, tf, coeffs, n_h, VR, p)
    ax = coeffs(:,1); ay = coeffs(:,2); az = coeffs(:,3);

    x = ax(1); y = ay(1); z = az(1);
    xdot = 0; ydot = 0; zdot = 0; xddot = 0; yddot = 0; zddot = 0;
    for i = 2:n_h+1
        x = x + ax(i)*cos(2*pi*(i-1)*t/tf) + ax(i+n_h)*sin(2*pi*(i-1)*t/tf);
        y = y + ay(i)*cos(2*pi*(i-1)*t/tf) + ay(i+n_h)*sin(2*pi*(i-1)*t/tf);
        z = z + az(i)*cos(2*pi*(i-1)*t/tf) + az(i+n_h)*sin(2*pi*(i-1)*t/tf);
        xdot = xdot - ax(i)*(2*pi*(i-1)/tf)*sin(2*pi*(i-1)*t/tf) + ax(i+n_h)*(2*pi*(i-1)/tf)*cos(2*pi*(i-1)*t/tf);
        ydot = ydot - ay(i)*(2*pi*(i-1)/tf)*sin(2*pi*(i-1)*t/tf) + ay(i+n_h)*(2*pi*(i-1)/tf)*cos(2*pi*(i-1)*t/tf);
        zdot = zdot - az(i)*(2*pi*(i-1)/tf)*sin(2*pi*(i-1)*t/tf) + az(i+n_h)*(2*pi*(i-1)/tf)*cos(2*pi*(i-1)*t/tf);
        xddot = xddot - ax(i)*((2*pi*(i-1)/tf)^2)*cos(2*pi*(i-1)*t/tf) - ax(i+n_h)*((2*pi*(i-1)/tf)^2)*sin(2*pi*(i-1)*t/tf);
        yddot = yddot - ay(i)*((2*pi*(i-1)/tf)^2)*cos(2*pi*(i-1)*t/tf) - ay(i+n_h)*((2*pi*(i-1)/tf)^2)*sin(2*pi*(i-1)*t/tf);
        zddot = zddot - az(i)*((2*pi*(i-1)/tf)^2)*cos(2*pi*(i-1)*t/tf) - az(i+n_h)*((2*pi*(i-1)/tf)^2)*sin(2*pi*(i-1)*t/tf);
    end
    
    sigma = [x,y,z,xdot,ydot,zdot,xddot,yddot,zddot];
    [state, control] = get_xu(sigma, p, VR);
    state = [state, [x, y, z]];
end

function [x,u] = get_xu(sigma, p, VR)
    m = 4.5;
    rho = 1.225;
    S = 0.473;
    g = 9.806;
    Cd0 = 0.0173;
    Cd1 = -0.0337;
    Cd2 = 0.0517;
    
    z = sigma(3);
    zdot = sigma(6); xdot = sigma(4); ydot = sigma(5); 
    zddot = sigma(9); xddot = sigma(7); yddot = sigma(8);

    % wind model
    p_exp = p;
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
    
    % aerodynamic forces
    Cd = Cd0+ Cd1*Cl + Cd2*Cl^2;
    D = 0.5*rho*S*V^2*Cd;
    T = m*Vdot + D + m*g*sin(gamma) + m*Wxz*zdot*cos(gamma)*cos(chi);
    CT = T/(0.5*rho*S*V^2);
    x = [V, chi, gamma]; u = [Cl, nu, CT];
end