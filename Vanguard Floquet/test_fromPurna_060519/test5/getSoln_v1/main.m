clear all

% 19/04/19
% van der pol oscillator
% xdot = - y
% ydot = - mu(x^2 - 1)*y - x

N = 20;

mu = 0.9;

x0 = struct; x0.val = 2; x0.swth = 0; 

Tinit = 6;
t = linspace(0,Tinit,N+1)'; t(1) = [];
Zinit = [x0.val*cos(2*pi*t/Tinit);-(2*pi*x0.val/Tinit)*sin(2*pi*t/Tinit);Tinit];


[D,~] = fourierDiff(N);

% generate guess from previous solution
clear('Zinit');
load results;
Zinit(2*N+1,1) = 0;
Nold = (length(Z)-1)/2;
Zinit(1:N) = sincInterp(linspace(0,1,Nold),Z(1:Nold),linspace(0,1,N));
Zinit(N+1:2*N) = sincInterp(linspace(0,1,Nold),Z(Nold+1:2*Nold),linspace(0,1,N));
Zinit(end) = Z(end);
clear('Z','Nold');

alg = 'levenberg-marquardt';
opts = optimoptions('lsqnonlin','Algorithm',alg,'Display','iter');
[Z,resnorm] = lsqnonlin(@(Z) residFunc(Z,D,N,mu,x0),Zinit,[],[],opts);

T = Z(end);
figure
plot([Z(N);Z(1:N)],[Z(2*N);Z(N+1:2*N)],'-r','LineWidth',1);
hold on
plot(Z(1),Z(N+1),'ob');
axis equal

save('results','Z');
