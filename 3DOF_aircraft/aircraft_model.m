function Ydot = aircraft_model(t, Y, u, par)
    Cl = u(1); nu = u(2);
    V = Y(1,1); psi = Y(2,1); gamma = Y(3,1); h = Y(4,1); x = Y(5,1); y = Y(6,1);
    % model_par = [m, rho, S, g, Cd0, Cd1, Cd2, b, hR, VR] 
    m = par(1); rho = par(2); S = par(3); g = par(4); Cd0 = par(5); Cd1 = par(6); Cd2 = par(7);
    b = par(8); hR = par(9); VR = par(10);
    
    p = 0.143;
    
    Cd = Cd0 + Cd1*Cl + Cd2*Cl^2;
    L = 0.5*rho*S*Cl*V^2;
    D = 0.5*rho*S*Cd*V^2;
    Vw = wind_model(h, [VR, hR, p,1]);
    
    hdot = V*sin(gamma);
    xdot = V*cos(gamma)*sin(psi) + Vw;
    ydot = V*cos(gamma)*cos(psi);
    
    Vwdot = VR*hdot;
    
    Vdot = -(D/m) - g*sin(gamma) - Vwdot*cos(gamma)*sin(psi);
    psidot = (L*sin(nu) - m*Vwdot*cos(psi))/(m*V*cos(gamma));
    gammadot = (L*cos(nu) - m*g*cos(gamma) + m*Vwdot*sin(gamma)*sin(psi))/(m*V);
    
    Ydot = [Vdot; psidot; gammadot; hdot; xdot; ydot];
end