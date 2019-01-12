clc

addpath('../trajectory')
addpath('../constraint_funcs')
addpath('../floquet')
addpath('../solutions')
addpath('../aircraft')

load('../solutions/lin_O.mat');

% ac = aircraft(); ac.tf = t(end); ac.N = 8; ac.VR = VR;
% ac.coeffs = get_coeffs([x,y,z], ac.tf, ac.N);
t = linspace(0, solution.tf, 10000);
state_DF = zeros(6, length(t));
for i = 1:length(t)
    sig = get_traj(t(i), ac.tf, ac.coeffs, ac.N);
    ac = ac.get_xu(sig);
    state_DF(:,i) = [ac.x(1);ac.x(3);ac.x(2);sig(1);sig(2);sig(3)];
end

sig_0 = get_traj(t(1), ac.tf, ac.coeffs, ac.N);
ac = ac.get_xu(sig_0);
y0 = [ac.x(1);ac.x(3);ac.x(2);sig(1);sig(2);sig(3)];
[t, state] = ode45(@(t,y) ac.non_flat_model(t, y), t, y0);

choice = 5;
plot(t, state_DF(choice,:), '--b', 'LineWidth', 1.5);
hold on; grid on;
plot(t, state(:,choice), '-r', 'LineWidth', 0.75);

error = state_DF - state';

rmpath('../aircraft')
rmpath('../trajectory')
rmpath('../constraint_funcs')
rmpath('../floquet')
rmpath('../solutions')
