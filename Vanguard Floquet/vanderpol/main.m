clearvars
d = 2;
global mu;
%% mu = - 0.1
load('mu_m01.mat')
[~,x] = fourierdiff(0.5*(length(X)-1));
tt = x*X(end,1)/(2*pi); % grid over which periodic solution was generated
N = 50;
lam = spectralFE(N, d, tt, X);
subplot(311)
scatter(real(lam), imag(lam), 'ob', 'DisplayName', strcat('N = ', num2str(N)));
hold on
N = 500;
lam = spectralFE(N, d, tt, X);
scatter(real(lam), imag(lam), 'xm', 'DisplayName', strcat('N = ', num2str(N)));
FM = retrieve_FM(d, tt, X);
scatter(real(FM), imag(FM), '*r', 'DisplayName', 'True FM');
theta = linspace(0, 2*pi, 1000);
% plot(cos(theta), sin(theta), 'k', 'LineWidth', 1.25, 'HandleVisibility', 'off');
title(strcat('Variation with N for $\mu = $', num2str(mu)), 'Interpreter', 'latex');
legend; grid minor
%% mu = + 0.1
load('mu_01.mat')
[~,x] = fourierdiff(0.5*(length(X)-1));
tt = x*X(end,1)/(2*pi); % grid over which periodic solution was generated
N = 50;
lam = spectralFE(N, d, tt, X);
subplot(312)
scatter(real(lam), imag(lam), 'ob', 'DisplayName', strcat('N = ', num2str(N)));
N = 500;
hold on
lam = spectralFE(N, d, tt, X);
scatter(real(lam), imag(lam), 'xm', 'DisplayName', strcat('N = ', num2str(N)));
FM = retrieve_FM(d, tt, X);
scatter(real(FM), imag(FM), '*r', 'DisplayName', 'True FM');
theta = linspace(0, 2*pi, 1000);
% plot(cos(theta), sin(theta), 'k', 'LineWidth', 1.25, 'HandleVisibility', 'off');
title(strcat('Variation with N for $\mu = $', num2str(mu)), 'Interpreter', 'latex');
legend; grid minor
%% mu = + 1
load('mu_1.mat')
[~,x] = fourierdiff(0.5*(length(X)-1));
tt = x*X(end,1)/(2*pi); % grid over which periodic solution was generated
N = 50;
lam = spectralFE(N, d, tt, X);
subplot(313)
scatter(real(lam), imag(lam), 'ob', 'DisplayName', strcat('N = ', num2str(N)));
N = 500;
hold on
lam = spectralFE(N, d, tt, X);
scatter(real(lam), imag(lam), 'xm', 'DisplayName', strcat('N = ', num2str(N)));
FM = retrieve_FM(d, tt, X);
scatter(real(FM), imag(FM), '*r', 'DisplayName', 'True FM');
theta = linspace(0, 2*pi, 1000);
% plot(cos(theta), sin(theta), 'k', 'LineWidth', 1.25, 'HandleVisibility', 'off');
title(strcat('Variation with N for $\mu = $', num2str(mu)), 'Interpreter', 'latex');
legend;
grid minor
%% mu = - 4
% load('mu_m1.mat')
% [~,x] = fourierdiff(0.5*(length(X)-1));
% tt = x*X(end,1)/(2*pi); % grid over which periodic solution was generated
% % N = 50;
% % lam = spectralFE(N, d, tt, X);
% % subplot(224)
% % scatter(real(lam), imag(lam), 'ob', 'DisplayName', strcat('N = ', num2str(N)));
% N = 500;
% hold on
% lam = spectralFE(N, d, tt, X);
% scatter(real(lam), imag(lam), 'xm', 'DisplayName', strcat('N = ', num2str(N)));
% FM = retrieve_FM(d, tt, X);
% scatter(real(FM), imag(FM), '*r', 'DisplayName', 'True FM');
% theta = linspace(0, 2*pi, 1000);
% % plot(cos(theta), sin(theta), 'k', 'LineWidth', 1.25, 'HandleVisibility', 'off');
% title(strcat('Variation with N for $\mu = $', num2str(mu)), 'Interpreter', 'latex');
% legend; axis equal;
%% functions 
function lam = spectralFE(N, d, tt, X)
    [D,x] = fourierdiff(N);
    T = X(end,1);
    t = x/(2*pi/T);
    M = zeros(N*d, N*d);
    DD = zeros(N*d, N*d);
    for i = 1:length(t)
        A = vanderpol_jac(t(i), tt, X);
        for k = 1:d
            kk = (k-1)*N;
            for j = 1:d
                jj = (j-1)*N;
                M(i+kk,i+jj) = A(k,j);
            end
        end
    end
    for i = 1:d
        j = (i-1)*N+1;
        DD(j:j+N-1,j:j+N-1) =2*pi*D/T;
    end
    lam = -eig(DD-M);
%     lam = exp(lam*T);
end
function A = vanderpol_jac(t, tt, X)
    global mu;
    N = 0.5*(length(X)-1);
    x = interp_sinc(tt, X(1:N,1), t);
    xdot = interp_sinc(tt, X(N+1:2*N,1), t);
    A = zeros(2);
    A(1,2) = 1;
    A(2,1) = -1 - 2*mu*x*xdot; A(2,2) = mu*(1-x^2);
end
function FM = retrieve_FM(d, tt, X)
    T = X(end,1);
    % FTM by explicit time evolution (N-pass)
    y_0 = eye(d); tspan = [0, T]; y_T = zeros(d);
    options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);
    for i = 1:d
        [~,y] = ode45(@(t,y) linear_vdp(t, y, tt, X), tspan, y_0(:,i), options);
        y_T(:,i) = y(end,:)';
    end
    FTM = y_T;
    FM = eig(FTM);
    FM = log(FM)/T;
end
function dy = linear_vdp(t, y, tt, X)
    A = vanderpol_jac(t, tt, X);
    dy = A*y;
end