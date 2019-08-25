function Ydot = dynamic_model_DS(x, u, p, VR)
    m = 4.5;
    rho = 1.225;
    S = 0.473;
    g = 9.806;
    Cd0 = 0.0173;
    Cd1 = -0.0337;
    Cd2 = 0.0517;
    
    V = x(1,1); gamma = x(3,1); chi = x(2,1);
    z = x(6,1);  
    Cl = u(1,1); nu = u(2,1); CT = 0;

    Cd = Cd0 + Cd1*Cl + Cd2*Cl^2;
    D = 0.5*rho*S*V^2*Cd;
    L = 0.5*rho*S*V^2*Cl;
    T = 0.5*rho*S*V^2*CT;

    % wind model
    Wx = VR*(-z)^p;
    Wxz = (p*VR)*((-z)^p)/z;

    zdot = -V*sin(gamma); xdot = V*cos(chi)*cos(gamma) + Wx;
    ydot = V*sin(chi)*cos(gamma);
    Vdot = (-D/m) - g*sin(gamma) - Wxz*zdot*cos(chi)*cos(gamma) + (T/m);
    chidot = (L*sin(nu)/(m*V*cos(gamma))) + Wxz*zdot*sin(chi)/(V*cos(gamma));   
    gammadot = (L*cos(nu)/(m*V)) - g*cos(gamma)/V + Wxz*zdot*cos(chi)*sin(gamma)/V;

    Ydot = [Vdot; chidot; gammadot; xdot; ydot; zdot];
end