N = 8;

aircraft = aircraft();
aircraft.limits
[coeffs, VR, tf, sol] = generate_traj(aircraft, N, []);
plot_traj(coeffs, tf, N);
%[coeffs, VR, tf, xguess] = generate_traj(aircraft, N, sol);

%x_guess = x;

%[a, eta, tf, VR, eig_vec, eig_val, x] = get_floquet(limits, model_par, x_guess, N, M);