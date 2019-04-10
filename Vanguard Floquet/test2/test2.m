% test hypothesis on Chebfun 
% http://www.chebfun.org/examples/ode-linear/Floquet.html

clear all

global a d
a = 0.15;
d = 4;
N = 4;

[D, x] = fourierDiff(N);
T = pi; % pi - periodic system
t = x*T/(2*pi);

Dmat = zeros(N*d);
M = zeros(d*N);
for i = 1:N
    A = sysModel(t(i));
    for j = 1:d
       for k = 1:d
            M(i+(j-1)*N,(k-1)*N+i) = A(j,k); 
       end
       Dmat((j-1)*N+1:j*N,(j-1)*N+1:j*N) = D*(2*pi/T);
    end
end

[eigVec,eigE] = eig(Dmat - M);
eigE = -1*diag(eigE);
eigM = exp(eigE*T);
FTM = get_FTM(t,0);
FM = eig(FTM);
FE = log(FM)/T;

%%

figure
hold on
thet = linspace(0,2*pi,100);
plot(cos(thet),sin(thet),'-k','LineWidth',1);
axis equal
idx = cell(d);
ix = [];
for i = 1:d
    idx{i} = find(abs(FM(i)-eigM)<0.05);
    ix = [ix;idx{i}];
end
ix = unique(ix);
ixx = 1:length(eigM);
ixx(ix) = [];  
p1 = plot(eigM(ix),'ob','LineWidth',1);
p0 = plot(eigM(ixx),'sm','LineWidth',1);
p2 = plot(real(FM),imag(FM),'*r','LineWidth',1);
legend([p1,p0,p2],{horzcat('Accurate FM from spectral ','(',num2str(length(ix)),')'),horzcat('Spurious FM from spectral ','(',num2str(length(ixx)),')'),'FM from time-evolve'});
title({'Coupled Oscillator',horzcat('N = ',num2str(N))});

function A = sysModel(t)
global a
    A = [             0, 1,               0, 0;
        -(2+a*cos(2*t)), 0,               1, 0;
                      0, 0,               0, 1;
                      1, 0, -(2+a*cos(2*t)), 0];
end

function FTM = get_FTM(t,flg)
global d
    switch flg
        case 0 % time-evolve method
            Phi=zeros(d);
            Phi0=eye(d);
            options = odeset('RelTol',1e-9,'AbsTol',1e-9);
            linDervs = @(t,X) sysModel(t)*X;
            for i=1:d
                [~,X]=ode15s(@(t,X) linDervs(t,X),[0,t(end)],Phi0(:,i),options);
                Phi(:,i)=X(end,:)';
            end
            FTM = Phi;
        case 1 % Friedmann's method
            h = t(2) - t(1);
            FTM = eye(d);
            for i = 1:length(t)-2
                K = friedmann_K(h, (t(end) - i*h));
                FTM = FTM*K;
            end
    end
 
end

% to calculate FTM
function K = friedmann_K(h, t)
global d
    A_psi = sysModel(t);
    A_psi_h_by_2 = sysModel(t+0.5*h);
    A_psi_h = sysModel(t+h);
    E = A_psi_h_by_2*(eye(d) + 0.5*h*A_psi);
    F = A_psi_h_by_2*(eye(d) + (-0.5 + 2^(-0.5))*h*A_psi + (1 - 2^(-0.5))*h*E);
    G = A_psi_h*(eye(d) - h*(2^(-0.5))*E + (1 + 2^(-0.5))*h*F);       
    K = eye(d) + (h/6)*(A_psi + 2*(1 - 2^(-0.5))*E + 2*(1 + 2^(-0.5))*F + G);
end