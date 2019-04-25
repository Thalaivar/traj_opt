clearvars
% load periodic solution
load('vdp1.mat')
N = 0.5*(length(X)-1);
T = X(2*N+1,1);
[~,x] = fourierdiff(N);
tt = x*T/(2*pi);
x1 = X(1:N,1); x2 = X(N+1:2*N,1);
t = linspace(0, T, 1000);
x1_interp = interp_sinc(tt, x1, t);
x2_interp = interp_sinc(tt, x2, t);

subplot(211)
plot(t, x1_interp, '-r', 'LineWidth', 1.25, 'DisplayName', 'interpolated');
hold on
plot(tt, x1, '-ob', 'LineWidth', 1.0, 'DisplayName', 'true');
xlabel('t'); ylabel('$x_1$', 'Interpreter', 'latex');
subplot(212)
plot(t, x2_interp, '-r', 'LineWidth', 1.25, 'DisplayName', 'interpolated');
hold on
plot(tt, x2, '--b', 'LineWidth', 1.0, 'DisplayName', 'true');
xlabel('t'); ylabel('$x_2$', 'Interpreter', 'latex');
legend;