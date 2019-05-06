clear all

% 11/04/19
% glycolysis oscillator
% xdot = -x + ay + x^2y
% ydot = b - ay - x^2y

N = 140;

Zinit = [0.5*ones(2*N,1);14];
a = 0.04;
b = 0.4;
x0 = 1;
[D,~] = fourierDiff(N);

% % generate guess from previous solution
% clear('Zinit');
% load results;
% Zinit(2*N+1,1) = 0;
% Nold = (length(Z)-1)/2;
% Zinit(1:N) = sincInterp(linspace(0,1,Nold),Z(1:Nold),linspace(0,1,N));
% Zinit(N+1:2*N) = sincInterp(linspace(0,1,Nold),Z(Nold+1:2*Nold),linspace(0,1,N));
% Zinit(end) = Z(end);
% clear('Z','Nold');

alg = 'levenberg-marquardt';
opts = optimoptions('lsqnonlin','Algorithm',alg,'Display','iter');
[Z,resnorm] = lsqnonlin(@(Z) residFunc(Z,D,N,a,b,x0),Zinit,[],[],opts);

T = Z(end);
figure
plot([Z(N);Z(1:N)],[Z(2*N);Z(N+1:2*N)],'-r','LineWidth',1);
hold on
plot(Z(1),Z(N+1),'ob');

save('results','Z');
