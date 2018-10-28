% get the trajectory info
clear all
load('/Users/dhruvlaad/ranjit_mohan/traj_opt/3DOF_aircraft/solutions/lin_O_shaped.mat')

aircraft = aircraft();
aircraft.traj_params.tf = tf; aircraft.traj_params.VR = VR;
aircraft.traj_params.coeffs = coeffs;

y_0 = eye(6); tspan = [0, tf]; y_T = zeros(6,6);
for i = 1:6
    [~,y] = ode45(@(t,y) model(t, y, aircraft, N), tspan, y_0(:,i));
    y_T(:,i) = y(end,:)';
end

H = y_T;

function ydot = model(t, y, aircraft, N)
    tf = aircraft.traj_params.tf; 
    VR = aircraft.traj_params.VR;
    coeffs = aircraft.traj_params.coeffs;
    
    A = aircraft.get_jac(t, tf, VR, coeffs, N);
    ydot = A*y;
end