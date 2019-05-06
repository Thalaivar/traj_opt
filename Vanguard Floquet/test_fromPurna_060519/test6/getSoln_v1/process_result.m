clear all


% load('./results/result_A-0p04-0p4')
% load('./results/result_A-0p02-0p8')
% load('./results/result_A-0p02-0p2')
% load('./results/result_A-0p11-0p6')
% load('./results/result_A-0p08-0p8')
load('./results/result_A-0p01-0p8')

T = Z(end);
N = (length(Z)-1)/2;

solnData = struct;
solnData.x1 = Z(1:N);
solnData.x2 = Z(N+1:2*N);
solnData.A = A;
solnData.T = T;
solnData.t = linspace(0,T,N+1); solnData.t(1) = [];

tt = linspace(0,T,1000); tt(1) = [];
xx1 = sincInterp(solnData.t,solnData.x1,tt);
xx2 = sincInterp(solnData.t,solnData.x2,tt);

figure
plot(solnData.x1,solnData.x2,'-r');
hold on
plot(xx1,xx2,'--b');
titl_str = horzcat('Phase Plane (a = ',num2str(solnData.A(1)),', b = ',num2str(solnData.A(2)),')');
title(titl_str);
axis equal
grid minor

figure
plot(solnData.t,solnData.x1,'-b');

save('./results/solnData_A-0p01-0p8','solnData');