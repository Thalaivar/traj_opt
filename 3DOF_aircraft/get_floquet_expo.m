load('/Users/dhruvlaad/ranjit_mohan/traj_opt/3DOF_aircraft/solutions/lin_O_shaped.mat')
M = 100;
aircraft = aircraft();
aircraft.traj_params.tf = tf; aircraft.traj_params.VR = VR;
aircraft.traj_params.coeffs = coeffs;
[eig_vec, eig_val, solution] = get_floquet(aircraft, N, M, []);