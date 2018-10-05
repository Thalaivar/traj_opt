function [V, gamma, psi, nu, Cl, T, CT] = DF_aircraft_model(z, wind_par, model_par)

    h = z(1);
    hdot = z(4); xdot = z(5); ydot = z(6); 
    hddot = z(7); xddot = z(8); yddot = z(9);
    VR = wind_par(1); choice = wind_par(2);
    m = model_par(1); rho = model_par(2); S = model_par(3); g = model_par(4);
    Cd0 = model_par(5); Cd1 = model_par(6); Cd2 = model_par(7); b = model_par(8);
    
    a = 0.143; hR = 20; % wind model param
    Vw = wind_model(h, [VR, hR, a, choice]);
    
    if choice == 1
        Vwdot = VR*hdot;
    else
        Vwdot = VR*a*Vw*hdot/h;
    end
    
    V        = (hdot^2 + (xdot - Vw)^2 + ydot^2)^0.5;
    Vdot     = (hdot*hddot + (xdot - Vw)*(xddot - Vwdot) + ydot*yddot)/V;
    gamma    = asin(hdot/V);
    gammadot = (V*hddot - hdot*Vdot)/(V*(V^2 - hdot^2)^0.5);
    psi      = atan((xdot - Vw)/ydot);
    psidot   = (ydot*(xddot - Vwdot) - (xdot - Vw)*yddot)/(ydot^2 + (xdot - Vw)^2);
    nu       = atan((cos(gamma)*psidot + (Vwdot/(V*cos(psi))))/(gammadot + (g/(V*cos(gamma))) - (Vwdot/(V*sin(gamma)*sin(psi)))));
    Cl       = (m*V*gammadot + m*g*cos(gamma) - m*Vwdot*sin(gamma)*sin(psi))/(0.5*rho*S*(V^2)*cos(nu));
    Cd       = Cd0 + Cd1*Cl + Cd2*Cl^2;
    D        = 0.5*rho*S*Cd*V^2;
    T        = m*Vdot + D + m*g*sin(gamma) + m*Vwdot*cos(gamma)*sin(psi);
    CT = T/(0.5*rho*S*V^2);
end