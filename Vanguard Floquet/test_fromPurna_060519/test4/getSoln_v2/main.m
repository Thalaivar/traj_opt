clear all
% 11/04/19
% glycolysis oscillator
% xdot = -x + ay + x^2y
% ydot = b - ay - x^2y

% % #1
% Z0init = [0.5;0.5;14];
% a = 0.04;
% b = 0.4;

% % #2
Z0init = [1;0.5;14];
a = 0.02;
b = 0.5;

% % #3
% Z0init = [1;1;14];
% a = 0.04;
% b = 0.7;

opts = odeset('RelTol',1e-11,'AbsTol',1e-11);
options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Display','iter');
[Z,resnorm] = lsqnonlin(@(Z) residFunc(Z,a,b,opts),Z0init,[],[],options);

T = Z(end);
[t,z] = ode15s(@(t,z) dynFunc(t,z,a,b),[0,T],Z(1:2),opts);

%%

figure
plot(z(:,1),z(:,2),'-r','LineWidth',1);
title('x_1 - x_2');

figure
hold on
plot(t,z(:,1),'-r','LineWidth',1);
plot(t,z(:,2),'-b','LineWidth',1);

save('soln2','t','z','a','b','T');