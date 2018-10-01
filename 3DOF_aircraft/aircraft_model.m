function Ydot = aircraft_model(t, Y, u, par)
    Cl = u(1,1); nu = u(2,1); T = u(3,1);
    V = Y(1,1); gamma = Y(2,1); psi  = Y(3,1); h = Y(4,1); x = Y(5,1); y = Y(6,1);
    rho = par(1); Cd0 = par(2); Cd1 = par(3); Cd2 = par(4); S = par(5); m = par(6); g = par(7);
    VR = par(8); hR = par(9);
    
    p = 0.143; % parameter for wind model
    
    Cd = Cd0 + Cd1*Cl + Cd2*Cl^2;
    L = 0.5*rho*S*Cl*V^2; D = 0.5*rho*S*Cd*V^2;
    Vw = wind_model(h, [VR, hR, p]);
    
    hdot = V*sin(gamma);
    xdot = V*cos(gamma)*sin(psi) + Vw;
    ydot = V*cos(gamma)*cos(psi);
    
    Vwdot = p*Vw*hdot/h;
    Vdot = (T - D - m*g*sin(gamma) - m*Vwdot*cos(gamma)*sin(psi))/m;
    psidot = (L*sin(nu) - m*Vwdot*cos(psi))/(m*V*cos(gamma));
    gammadot = (L*cos(nu) - m*g*cos(gamma) + m*Vwdot*sin(gamma)*sin(psi))/(m*V);
    
    Ydot = [Vdot, gammadot, psidot, hdot, xdot, ydot];
end