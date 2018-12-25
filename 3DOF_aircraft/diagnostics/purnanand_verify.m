clearvars

addpath('trajectory')
addpath('constraint_funcs')
addpath('floquet')
addpath('solutions')

load('3DOF_eight.mat');
ac = aircraft(); ac.tf = t(end); ac.VR = VR; ac.N = 8;

% retrieve coeffs
ac.coeffs = get_coeffs([x,y,z], ac.tf, ac.N);

tspan = [ac.tf/4, 10*ac.tf];
sig_0 = get_traj(tspan(1), ac.tf, ac.coeffs, ac.N);
ac = ac.get_xu(sig_0);
y0 = [ac.x(1);ac.x(3);ac.x(2);sig_0(1);sig_0(2);sig_0(3)];

[t, y] = ode45(@(t,y) ac.non_flat_model(t, y), tspan, y0);

t2 = linspace(0, ac.tf, 1000);
nom_traj = zeros(3, length(t2));

for i = 1:length(t2)
    sig = get_traj(t2(i), ac.tf, ac.coeffs, ac.N);
    nom_traj(:,i) = [sig(1);-sig(2);-sig(3)];
end

plot3(nom_traj(1,:), nom_traj(2,:), nom_traj(3,:), '--b', 'LineWidth', 1.5);
hold on
comet3(y(:,4), -y(:,5), -y(:,6))
grid on

rmpath('trajectory')
rmpath('constraint_funcs')
rmpath('floquet')
rmpath('solutions')