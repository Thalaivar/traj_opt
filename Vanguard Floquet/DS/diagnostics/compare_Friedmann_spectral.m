addpath('..\lib\')

clearvars
close all
load('../3DOF/full_state_discrete/lineSep2.mat')
trajData.type = 'full-state';
trajData.shape = 'circle';

trajData = stateControlMat(sol, N, p, trajData);
%% spectral method
[spectralFE,  eigE, AM, groupSizes] = spectralMethod(trajData);
% [FE,groupSize,Dmat,Mmat,AM] = spectralMethodMax(trajData);
% %% Freidmann method
% [~,freidMannGrid] = fourierdiff(500);
% [~,freidmannFE] = freidmannMethod(trajData, freidMannGrid);
% %% time march method
% [FTM,timeMarchFE] = timeMarchMethod(trajData);
% %% compare difference in the two
% freidmannEstimate = sort(real(freidmannFE)', 'descend');
% timeMarchEstimate = sort(real(timeMarchFE)', 'descend');
% residualFreidmann = max(freidmannEstimate - timeMarchEstimate)
% %% plotting
% figure
% hold on
% p1 = plot(real(eigE),imag(eigE),'xm');
% p3 = plot(real(timeMarchFE), imag(timeMarchFE), 'xg', 'LineWidth', 1.25, 'MarkerSize', 10);
% p2 = plot(real(freidmannFE),imag(freidmannFE),'or','LineWidth',1,'MarkerSize',10);
% limSet = 1.2*max(abs(real(eigE)));
% xlim([-limSet,limSet]);
% if p == 1
%       traj_type = 'linear';
% else
%     traj_type = 'exponential';
% end
% title_str1 = ['Dynamic soaring in ', traj_type, ' profile. '];
% title_str2 = ['N = ' num2str(trajData.N)];
% title([title_str1 title_str2]);
% grid minor
% legend([p1,p2,p3],{'Spectral FE','Friedmann-FE','Time-evolved FE'});
% 
% scatter(real(spectralFE), imag(spectralFE), 75, 'sb', 'LineWidth', 1.25, 'DisplayName', 'Estimated FE')

%% functions
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