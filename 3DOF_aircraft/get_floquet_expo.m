load('/Users/dhruvlaad/ranjit_mohan/traj_opt/3DOF_aircraft/solutions/lin_O_shaped.mat')
M = 100 ;

 get_FTM
 eigvals = (1/tf)*log(eig(H2));
 eigval = [0.0174, 0.0204];
 init_guess_floq
 
 %[c ,ceq] = constFun([init_guess;eigval(1);eigval(2)], aircraft, N, M, true)

aircraft = aircraft();
aircraft.traj_params.tf = tf; aircraft.traj_params.VR = VR;
aircraft.traj_params.coeffs = coeffs;
[eig_vec, eig_val, solution] = get_floquet(aircraft, N, M, [init_guess;eigval(1);eigval(2)]);