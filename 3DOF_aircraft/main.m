load('params.mat')

N = 5;

limits = [Clmax, Vmax, nu_min, nu_max, CTmin, CTmax, hmin];
model_par = [m, rho, S, g, Cd0, Cd1, Cd2, b];

[a, eta, tf, VR, xguess] = generate_traj(limits, model_par, N, []);

plot_traj(a, eta, tf);  