clear all

% 30/04/19
% glycolytic oscillator
% xdot = -x + ay + x^2y
% ydot = b - ay - x^2y

N = 750;

% A = [a,b]
A = [0.01,0.8];

x0 = struct; x0.val = 2; x0.swth = 0; 

Tinit = 6;
t = linspace(0,Tinit,N+1)'; t(1) = [];
Zinit = [x0.val*cos(2*pi*t/Tinit);-(2*pi*x0.val/Tinit)*sin(2*pi*t/Tinit);Tinit];


[D,~] = fourierDiff(N);

% generate guess from previous solution
clear('Zinit');
load results;
% load('./results/result_A-0p04-0p4','Z')
Zinit(2*N+1,1) = 0;
Nold = (length(Z)-1)/2;
Zinit(1:N) = sincInterp(linspace(0,1,Nold),Z(1:Nold),linspace(0,1,N));
Zinit(N+1:2*N) = sincInterp(linspace(0,1,Nold),Z(Nold+1:2*Nold),linspace(0,1,N));
Zinit(end) = Z(end);
clear('Z','Nold');

alg = 'levenberg-marquardt';
opts = optimoptions('lsqnonlin','Algorithm',alg,'Display','iter');
[Z,resnorm] = lsqnonlin(@(Z) residFunc(Z,D,N,A,x0),Zinit,[],[],opts);

T = Z(end);
figure
plot([Z(N);Z(1:N)],[Z(2*N);Z(N+1:2*N)],'-r','LineWidth',1);
hold on
plot(Z(1),Z(N+1),'ob');
axis equal

save('results','Z');
