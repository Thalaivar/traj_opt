function plotDiffFlatSol(sol, p)
    N = (length(sol)-2)/3;
    x = zeros(N, 6); u = zeros(N,3);
    T = sol(3*N+1); VR = sol(3*N+2);
    x(:,4:6) = [sol(1:N), sol(N+1:2*N), sol(2*N+1:3*N)];
    
    if p == 1
        profile = 'linear';
    else
        profile = 'exponential';
    end
    
    [D, t] = fourierdiff(N);
    t = T*t/(2*pi);
    diffFac = (2*pi/T);
    
    column2 = [-(N^2)/12-1/6, -((-1).^(1:(N-1)))./(2*(sin((1:(N-1))*pi/N)).^2)];
    DD = toeplitz(column2,column2([1, N:-1:2])); % second derivative matrix
    dx  = diffFac*D*x(:,4:6);
    ddx = diffFac*diffFac*DD*x(:,4:6);
    
    for i = 1:N
        [state, control] = diffFlatModel([x(i,4:6);dx(i,:);ddx(i,:)], p, VR);
        x(i,1:3) = state; u(i,:) = control;
    end
    
    subplot(311)
    plot(t, x(:,1)); xlabel('t (s)'); ylabel('$V$', 'Interpreter', 'latex');
    grid minor; title(['State time history for N = ', num2str(N), ' in ', profile, ' wind'])
    subplot(312)
    plot(t, x(:,2)); xlabel('t (s)'); ylabel('$\chi$', 'Interpreter', 'latex');
    grid minor;
    subplot(313)
    plot(t, x(:,3)); xlabel('t (s)'); ylabel('$\gamma$', 'Interpreter', 'latex');
    grid minor;
    figure
    plot3(x(:,4), -x(:,5), -x(:,6), '--om'); xlabel('x'); ylabel('y'); zlabel('z');
    title(['Trajectory for N = ', num2str(N), ' in ', profile, ' wind'])
    grid minor;
    figure
    subplot(311)
    plot(t, u(:,1)); xlabel('t (s)'); ylabel('$C_L$', 'Interpreter', 'latex');
    grid minor; title(['Control time history for N = ', num2str(N), ' in ', profile, ' wind'])
    subplot(312)
    plot(t, u(:,2)); xlabel('t (s)'); ylabel('$\mu$', 'Interpreter', 'latex');
    grid minor;
    subplot(313)
    plot(t, u(:,3)); xlabel('t (s)'); ylabel('$C_T$', 'Interpreter', 'latex');
    grid minor;
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