addpath('trajectory')
addpath('constraint_funcs')
addpath('floquet')
addpath('solutions')

load('lin_eight.mat')

tspan = [0, 10*tf];

% get control history
t = linspace(tspan(1), tspan(2), 1000);
U = zeros(3, length(t));
for i = 1:length(t)
    sig = get_traj(t(i), ac.tf, ac.coeffs, ac.N);
    ac = ac.get_xu(sig);
    U(:,i) = [ac.u(1);ac.u(2);ac.u(3)];
end

% fit polynomial to control history
n_poly = 6; poly_p = zeros(3, n_poly+1);
for i = 1:3
    poly_p(i,:) = polyfit(t, U(i,:), n_poly);
end

rmpath('solutions')
rmpath('trajectory')
rmpath('constraint_funcs')
rmpath('floquet')