clearvars
% to plot both stable and unstable solutions

addpath('../trajectory');
addpath('../floquet');
addpath('..');

% stable solution
load('../solutions/stability_opt/lin_O_stable.mat');
tspan = [0, 10*ac.tf];
sig_0 = get_traj(tspan(1), ac.tf, ac.coeffs, ac.N);
ac = ac.get_xu(sig_0);
y0 = [ac.x(1); ac.x(3); ac.x(2); sig_0(1); sig_0(2); sig_0(3)];
[t, y] = ode45(@(t, y) ac.non_flat_model(t, y), tspan, y0);
X_stable = [y(:,4), -y(:,5), -y(:,6)];

% plot unstable solution
load('../solutions/stability_opt/lin_O_unstable.mat');
tspan = [0, 10*ac.tf];
sig_0 = get_traj(tspan(1), ac.tf, ac.coeffs, ac.N);
ac = ac.get_xu(sig_0);
y0 = [ac.x(1); ac.x(3); ac.x(2); sig_0(1); sig_0(2); sig_0(3)];
[t, y] = ode45(@(t, y) ac.non_flat_model(t, y), tspan, y0);
X_unstable = [y(:,4), -y(:,5), -y(:,6)];


% translate start of unstable to start of stable
trans = zeros(1,3);
for i = 1:3
    trans(1,i) = X_stable(1,i) - X_unstable(1,i);
end
for i = 1:3
    X_unstable(:,i) = X_unstable(:,i) + trans(1,i);
end

% wind profile
%x_wind = linspace(-120,100,100);

plot3(X_stable(:,1), X_stable(:,2), X_stable(:,3), 'b', X_unstable(:,1), X_unstable(:,2), X_unstable(:,3), 'r');
hold on
plot3(X_stable(1,1), X_stable(1,2), X_stable(1,3), 'sg');
xlabel('x'); ylabel('y'); zlabel('z'); 
grid on
axis equal

rmpath('..');
rmpath('../trajectory');
rmpath('../floquet');