function [c, ceq] = constraintFunctionDF(X, p, D, DD)
    N = (length(X)-2)/3;
    x = zeros(N, 6); u = zeros(N,3);
    T = X(3*N+1); VR = X(3*N+2);
    x(:,4:6) = [X(1:N), X(N+1:2*N), X(2*N+1:3*N)];
    
    diffFac = (2*pi/T);
    dx  = diffFac*D*x(:,4:6);
    ddx = diffFac*diffFac*DD*x(:,4:6);
    
    for i = 1:N
        [state, control] = diffFlatModel([x(i,4:6);dx(i,:);ddx(i,:)], p, VR);
        x(i,1:3) = state; u(i,:) = control;
    end
    
    % constraints on non-flat states
    % z constraints\
    hmin = 0.15; b = 3;
    c(1:N,1) = x(:,6) + 0.5*b*sin(u(:,2));
    % V constraints
    c(1*N+1:2*N,1) = 10 - x(:,1);
    c(2*N+1:3*N,1) = x(:,1) - 80;
    % gamma constraints
    c(3*N+1:4*N,1) = -pi/4 - x(:,3);
    c(4*N+1:5*N,1) = x(:,3) - pi/4;
    % CL constraints
    c(5*N+1:6*N,1) = -0.2 - u(:,1);
    c(6*N+1:7*N,1) = u(:,1) - 1.17;
    % nu constraints
    c(7*N+1:8*N,1) = -pi/3 - u(:,2);
    c(8*N+1:9*N,1) = u(:,2) - pi/3;
    % CT constraints
    c(9*N+1:10*N,1) = -1e-6 - u(:,3);
    c(10*N+1:11*N,1) = u(:,3) - 1e-6;
    % mu rate^2 constraint
    c(11*N+1:12*N,1) = diffFac*diffFac*DD*u(:,2) - 20*pi/180;
    c(12*N+1:13*N,1) = -20*pi/180 - diffFac*diffFac*DD*u(:,2);
    % CL rate constraint
    c(13*N+1:14*N,1) = diffFac*D*u(:,1) - 0.4;
    c(14*N+1:15*N,1) = -0.4 - diffFac*D*u(:,1);
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