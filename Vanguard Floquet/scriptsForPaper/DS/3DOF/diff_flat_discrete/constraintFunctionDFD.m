function [c, ceq] = constraintFunctionDFD(X, trajParams, windShear)
    X = X';
    if strcmp(windShear, 'same')
        N = (length(X)-1)/3;
        VR = trajParams.VR;
    else
        N = (length(X)-2)/3;
        VR = X(3*N+2);
    end
    T = X(3*N+1);
    
    diffFac = 2*pi/T;
    D = trajParams.D; DD = trajParams.DD;
    x = [X(1:N), X(N+1:2*N), X(2*N+1:3*N)];
    dx = [diffFac*D*x(:,1), diffFac*D*x(:,2), diffFac*D*x(:,3)];
    ddx = [(diffFac^2)*DD*x(:,1), (diffFac^2)*DD*x(:,2), (diffFac^2)*DD*x(:,3)];
    
    p = trajParams.p;
    state = zeros(N,3); U = zeros(N,3);
    for i = 1:N
        [xx, uu] = diffFlatModel([x(i,:);dx(i,:);ddx(i,:)], p, VR);
        state(i,:) = xx; U(i,:) = uu;
    end
    
    hmin = 0; b = 3;
    c(1:N,1) = x(:,3) + 0.5*b*abs(sin(U(:,2))) + hmin;
    c(end+1:end+N,1) = 10 - state(:,1); 
    c(end+1:end+N,1) = state(:,1) - 80;
    c(end+1:end+N,1) = -pi/4 - state(:,3);
    c(end+1:end+N,1) = state(:,3) - pi/4;
    c(end+1:end+N,1) = -0.2 - U(:,1);
    c(end+1:end+N,1) = U(:,1) - 1.17;
    c(end+1:end+N,1) = -pi/3 - U(:,2);
    c(end+1:end+N,1) = U(:,2) - pi/3;
    c(end+1:end+N,1) = -1e-4 - U(:,3);
    c(end+1:end+N,1) = U(:,3) - 1e-4;
    
    c(end+1:end+N) = -20*pi/180 - (diffFac^2)*DD*U(:,2);
    c(end+1:end+N) = (diffFac^2)*DD*U(:,2) - 20*pi/180;
    c(end+1:end+N) = -0.7 - (diffFac^2)*DD*U(:,1);
    c(end+1:end+N) = (diffFac^2)*DD*U(:,1) - 0.7;
    
    ceq = [];
end

function [x, u] = diffFlatModel(X, p, VR)
    g = 9.81;
    m = 4.5;
    rho = 1.225;
    S = 0.473;
    g = 9.806;
    Cd0 = 0.0173;
    Cd1 = -0.0337;
    Cd2 = 0.0517;
    
    z = X(1,3);
    zdot = X(2,3); xdot = X(2,1); ydot = X(2,2); 
    zddot = X(3,3); xddot = X(3,1); yddot = X(3,2);

    % wind model
    Wx = VR*(-z)^p;
    Wxz = (p*VR)*((-z)^p)/z;

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
    Cd = Cd0 + Cd1*Cl + Cd2*Cl^2;
    D = 0.5*rho*S*V^2*Cd;
    T = m*Vdot + D + m*g*sin(gamma) + m*Wxz*zdot*cos(gamma)*cos(chi);
    CT = T/(0.5*rho*S*V^2);
    
    x = [V, chi, gamma]; u = [Cl, nu, CT];
end