% test hypothesis on Chebfun 
% http://www.chebfun.org/examples/ode-linear/Floquet.html

clear all
addpath('../')

global d
d = 2;
N = 100;

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
FTM = get_FTM(t);
FE = log(eig(FTM))/T;
FM = eig(FTM);

figure
hold on
p1 = plot(real(eigE),imag(eigE),'xm');
p2 = plot(real(FE),imag(FE),'or','LineWidth',1,'MarkerSize',10);
limSet = 1.2*max(abs(real(eigE)));
xlim([-limSet,limSet]);
title_str1 = horzcat('Coupled Oscillators Floquent Exponents');
title_str2 = horzcat('N = ',num2str(N));
title({title_str1,title_str2});
grid minor
legend([p1,p2],{'Spectral FE','Time-evolved FE'});

rmpath('../')

function A = sysModel(t)
     A = [-1 + 1.5*(cos(t))^2, 1 - 1.5*sin(t)*cos(t);
         -1 - 1.5*sin(t)*cos(t), -1 + 1.5*(sin(t))^2];
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

function FTM = get_FTM(t)
global d
    h = t(2) - t(1);
    FTM = eye(d);
    for i = 1:length(t)-2
        K = friedmann_K(h, (t(end) - i*h));
        FTM = FTM*K;
    end
end