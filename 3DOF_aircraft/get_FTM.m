% get the trajectory info
load('/Users/dhruvlaad/ranjit_mohan/traj_opt/3DOF_aircraft/solutions/lin_eight.mat')

aircraft = aircraft();
aircraft.traj_params.tf = tf; aircraft.traj_params.VR = VR;
aircraft.traj_params.coeffs = coeffs;

y_0 = eye(3); tspan = linspace(0,tf,10000); y_T = zeros(3,3);
for i = 1:3
    [~,y] = ode45(@(t,y) model(t, y, aircraft, N), tspan, y_0(:,i));
    y_T(:,i) = y(end,:)';
end

H1 = y_T;

t = linspace(tf, 0, 1000);
Jk_0 = (aircraft.get_jac(t(1), tf, VR, coeffs, N) - aircraft.get_jac(t(2), tf, VR, coeffs, N))/2;
H2 = expm(Jk_0*(t(1)-t(2)));
for i = 2:length(t)-1
    Jk = (aircraft.get_jac(t(i), tf, VR, coeffs, N) + aircraft.get_jac(t(i+1), tf, VR, coeffs, N))/2;
    H2 = H2*expm(Jk*(t(i)-t(i+1)));
end


function ydot = model(t, y, aircraft, N)
    tf = aircraft.traj_params.tf; 
    VR = aircraft.traj_params.VR;
    coeffs = aircraft.traj_params.coeffs;
    
    A = aircraft.get_jac(t, tf, VR, coeffs, N);
    %A = [-1 + 1.5*(cos(t))^2, 1 - 1.5*sin(t)*cos(t); -1 - 1.5*sin(t)*cos(t), -1 + 1.5*(sin(t))^2];
    ydot = A*y;
end