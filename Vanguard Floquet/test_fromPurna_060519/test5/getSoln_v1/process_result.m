clear all

% load('./results/result_mu_p2_01')
% load('./results/result_mu_m2_01')

% load('./results/result_mu_p4_01')
% load('./results/result_mu_m4_01')

% load('./results/result_mu_m3_01')
% load('./results/result_mu_p3_01')

% load('./results/result_mu_p7_01')

% load('./results/result_mu_p0p9_M25')
% load('./results/result_mu_p0p9_M100')
% load('./results/result_mu_p0p9_M5')
load('./results/result_mu_p0p9_M10')

T = Z(end);
N = (length(Z)-1)/2;

solnData = struct;
solnData.x1 = Z(1:N);
solnData.x2 = Z(N+1:2*N);
solnData.mu = mu;
solnData.T = T;
solnData.t = linspace(0,T,N+1); solnData.t(1) = [];

tt = linspace(0,T,1000); tt(1) = [];
xx1 = sincInterp(solnData.t,solnData.x1,tt);
xx2 = sincInterp(solnData.t,solnData.x2,tt);

figure
plot(solnData.x1,solnData.x2,'-r');
hold on
plot(xx1,xx2,'--b');
titl_str = horzcat('Phase Plane (\mu = ',num2str(solnData.mu),')');
title(titl_str);
axis equal
grid minor

save('./results/solnData_mu_p0p9_M10','solnData');