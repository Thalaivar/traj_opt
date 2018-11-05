load('/Users/dhruvlaad/ranjit_mohan/traj_opt/3DOF_aircraft/solutions/lin_O_shaped.mat')

aircraft = aircraft();
aircraft.traj_params.tf = tf; aircraft.traj_params.VR = VR;
aircraft.traj_params.coeffs = coeffs;

[~, cheb_x] = cheb_diff(M-1);
cheb_t = 0.5*tf*(1 - cheb_x);
params = [eigval(1), eigval(2), N];

y0 = (1/sqrt(6))*ones(6,1);
tspan = cheb_t;

[t, y] = ode45(@(t,y) model(t, y, params, aircraft), tspan, y0);

init_guess = zeros(6*M,1);
for i = 1:M
    j = (i-1)*3 + 1;
    init_guess(j:j+2,1) = y(i,1:3);
    init_guess(3*M+j:3*M+j+2,1) = y(i,4:6);
end

function ydot = model(t, y, params, aircraft)
    alpha = y(1:3,1); beta = y(4:6,1);
    eigval = complex(params(1), params(2));
    N = params(3);
    tf = aircraft.traj_params.tf; VR = aircraft.traj_params.VR; coeffs = aircraft.traj_params.coeffs;
    A = aircraft.get_jac(t, tf, VR, coeffs, N);
    I = eye(3);
    alpha_dot = (A - real(eigval)*I)*alpha + imag(eigval)*beta;
    beta_dot  = (A - real(eigval)*I)*beta  - imag(eigval)*alpha; 
    ydot = [alpha_dot;beta_dot];
end