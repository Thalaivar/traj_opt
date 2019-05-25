clear all
close all
addpath('lib\')
load('full_state_discrete\solutions\maxObjFun\SEE60_sameWind.mat')
trajData.type = 'full-state'; trajData.shape = 'eight';
trajData = stateControlMat(sol, N, p, trajData);

[eigE, FE, maxViolation, grpSizes] = spectralMethod(trajData, false);
eigE = eigE((abs(imag(eigE)) <= 0.25*N*2*pi/trajData.T));
scatter(real(eigE), imag(eigE), 'xm')
% 
% [maxViolation, index] = max(FE(3,:));
% display(['Maximum relative violation in real part = ', num2str(maxViolation)])
% display(['Group ID for which maximum violation occurs : ', num2str(index)])
% FE
function trajData = stateControlMat(sol, N, p, trajData)
    x = zeros(N, 6);
    u = zeros(N, 3);
    if strcmp(trajData.type, 'full-state')
        for i = 1:8
            j = (i-1)*N + 1;
            if i <= 6
                x(:,i) = sol(j:j+N-1,1);
            else
                u(:,i-6) = sol(j:j+N-1,1);
            end
        end
        T = sol(8*N+1,1); VR = sol(8*N+2,1);
        [D,fourierGrid] = fourierdiff(N); t = T*fourierGrid/(2*pi);
        if strcmp(trajData.shape, 'circle')
            chiLinearTerm = sol(8*N+3);
            for i = 1:N
                x(i,2) = x(i,2) + chiLinearTerm*t(i);
            end
        end
    elseif strcmp(trajData.type, 'diff-flat')
        T = sol(3*N+1); VR = sol(3*N+2);
        diffFac = 2*pi/T;

        [D,fourierGrid] = fourierdiff(N);
        column2 = [-(N^2)/12-1/6, -((-1).^(1:(N-1)))./(2*(sin((1:(N-1))*pi/N)).^2)];
        DD = toeplitz(column2,column2([1, N:-1:2])); 
        x = [sol(1:N), sol(N+1:2*N), sol(2*N+1:3*N)];
        dx = [diffFac*D*x(:,1), diffFac*D*x(:,2), diffFac*D*x(:,3)];
        ddx = [(diffFac^2)*DD*x(:,1), (diffFac^2)*DD*x(:,2), (diffFac^2)*DD*x(:,3)];

        state = zeros(N,3); u = zeros(N,3);
        for i = 1:N
            [xx, uu] = diffFlatModel([x(i,:);dx(i,:);ddx(i,:)], p, VR);
            state(i,:) = xx; u(i,:) = uu;
        end
        x = [state, x];
    end
    trajData.T = T; trajData.VR = VR; trajData.fourierGrid = fourierGrid;
    trajData.X = x; trajData.U =  u; trajData.N = N; trajData.p = p; trajData.D = D;
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