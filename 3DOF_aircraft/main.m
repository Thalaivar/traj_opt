load('params.mat')

N = 5;
M = 50;
limits = [Clmax, Vmax, nu_min, nu_max, CTmin, CTmax, hmin];
model_par = [m, rho, S, g, Cd0, Cd1, Cd2, b];

%[~, ~, ~, ~, xguess] = generate_traj(limits, model_par, N, []);
%[a, eta, tf, VR, x] = generate_traj(limits, model_par, N, xguess);

%x_guess = x;

[a, eta, tf, VR, eig_vec, eig_val, x] = get_floquet(limits, model_par, x_guess, N, M);