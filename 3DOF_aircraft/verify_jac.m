clear all
load('/Users/dhruvlaad/ranjit_mohan/traj_opt/3DOF_aircraft/solutions/lin_O_shaped.mat')

J_cap = zeros(6);
t = linspace(0, tf, 10000);
aircarft = aircraft();
aircraft.traj_params.tf = tf; aircraft.traj_params.VR = VR;
aircraft.traj_params.coeffs = coeffs;
error = zeros(1,length(t)-1);
for k = 1:length(t)-1
    sigma = aircraft.get_traj(t(k), tf, coeffs, N);
    aircraft = aircraft.get_xu(sigma, VR);
    for i = 1:6
        x0 = [aircraft.x(1);aircraft.x(3);aircraft.x(2);sigma(1);sigma(2);sigma(3)];
        u = aircraft.u'; I = eye(6);
        h = t(k+1) - t(k);
        f1 = aircraft.model(x0-h*I(:,i),u, VR); f2 = aircraft.model(x0+h*I(:,i),u, VR);
        J_cap(:,i) = (f2 - f1)/(2*h);
    end
    
    A = aircraft.get_jac(t(k), tf, VR, coeffs, N);
    error(k) = norm(J_cap - A)/norm(A);
end
    