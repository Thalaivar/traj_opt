addpath('../trajectory')
addpath('../constraint_funcs')
addpath('../floquet')
addpath('../solutions')
addpath('..')

%load('../solutions/trajectory_opt/expo_O.mat')

% time evolution
tspan = [0,ac.tf];
sig_0 = get_traj(tspan(1), ac.tf, ac.coeffs, ac.N);
ac = ac.get_xu(sig_0);
y0 = [ac.x(1);ac.x(3);ac.x(2);sig_0(1);sig_0(2);sig_0(3)];
options = odeset('AbsTol', 1e-11, 'RelTol', 1e-11, 'Events', @(t,y) z_event(t,y));
sol = ode45(@(t,y) ac.non_flat_model(t, y), tspan, y0, options);

rmpath('../solutions')
rmpath('../trajectory')
rmpath('../constraint_funcs')
rmpath('../floquet')
rmpath('..')