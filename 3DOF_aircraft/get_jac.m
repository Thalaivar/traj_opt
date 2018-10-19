function A = get_jac(t, a, eta, tf, model_par, wind_par)
    ph = [a(:,1),[eta(:,1);0]]; px = [a(:,2),[eta(:,2);0]]; py = [a(:,3),[eta(:,3);0]];
    [h, x, y, hdot, xdot, ydot, hddot, xddot, yddot] = get_traj(t, ph, px, py, tf);
    z = [h, x, y, hdot, xdot, ydot, hddot, xddot, yddot];
    
    % state value at given point in trajectory
    [V, gamma, psi, nu, Cl, ~, ~] = DF_aircraft_model(z, wind_par, model_par);
    
    % model params
    m = model_par(1); rho = model_par(2); S = model_par(3); g = model_par(4);
    Cd0 = model_par(5); Cd1 = model_par(6); Cd2 = model_par(7);
    
    % wind value
    VR = wind_par(1); choice = wind_par(2);
    Vwdot = VR*hdot;
    
    % forces and drag polar
    Cd = Cd0 + Cd1*Cl + Cd2*Cl^2;
    L = 0.5*S*rho*Cl*V^2;
    
    % constructing the Jacobian 
    f11 = (1/m)*(-rho*S*Cd*V); f12 = (1/m)*(-m*Vwdot*cos(gamma)*cos(psi));
    f13 = (1/m)*(-m*g*cos(gamma) + m*Vwdot*sin(gamma)*sin(psi));
    f14 = 0; f15 = 0; f16 = 0;
    f21 = (1/(m*cos(gamma)))*((0.5*rho*S*Cl*sin(nu)) + m*Vwdot*cos(psi)/V^2);
    f22 = (1/(V*cos(gamma)))*(Vwdot*sin(psi));
    f23 = ((L*sin(nu) - m*Vwdot*cos(psi))/(m*V))*tan(gamma)*sec(gamma);
    f24 = 0; f25 = 0; f26 = 0; 
    f31 = (1/m)*(0.5*S*rho*Cl*cos(nu) + ((m*g*cos(gamma) - m*Vwdot*sin(gamma)*sin(psi))/V^2));
    f32 = Vwdot*sin(gamma)*cos(psi)/V;
    f33 = (g*sin(gamma) + Vwdot*cos(gamma)*sin(psi))/V;
    f34 = 0; f35 = 0; f36 = 0;
    f41 = sin(gamma); f42 = 0; f43 = V*cos(gamma); f44 = 0; f45 = 0; f46 = 0;
    f51 = cos(gamma)*sin(psi); f52 = V*cos(gamma)*cos(psi); f53 = -V*sin(gamma)*sin(psi);
    f54 = VR; f55 = 0; f56 = 0;
    f61 = cos(gamma)*cos(psi); f62 = -V*cos(gamma)*sin(psi); f63 = -V*sin(gamma)*cos(psi); 
    f64 = 0; f65 = 0; f66 = 0; 
    
    A = [f11, f12, f13, f14, f15, f16;
         f21, f22, f23, f24, f25, f26;
         f31, f32, f33, f34, f35, f36;
         f41, f42, f43, f44, f45, f46;
         f51, f52, f53, f54, f55, f56;
         f61, f62, f63, f64, f65, f66];
end