clear all
load('/Users/dhruvlaad/ranjit_mohan/traj_opt/3DOF_aircraft/solutions/lin_O_shaped.mat')

J_cap = zeros(3);
ac = aircraft();
ac.tf = tf; ac.VR = VR;
ac.coeffs = coeffs; ac.N = 8;

t = linspace(0, ac.tf, 1000);
error = zeros(1,length(t)-1);

for k = 1:length(t)-1
    sigma = get_traj(t(k), tf, coeffs, N);
    ac = ac.get_xu(sigma);
    for i = 1:3
        x0 = [ac.x(1);ac.x(3);ac.x(2)];
        u = ac.u'; I = eye(3);
        h = t(k+1) - t(k);
        f1 = ac.model(t(k), x0-h*I(:,i), u); f2 = ac.model(t(k), x0+h*I(:,i), u);
        J_cap(:,i) = (f2 - f1)/(2*h);
    end
    
    A = ac.get_jac(t(k));
    error(k) = norm(J_cap - A)/norm(A);
end

