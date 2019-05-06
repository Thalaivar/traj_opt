% 11/04/2019
clear all

addpath('./getSoln_v2');
load solnFourier3

global d
d = 2;
N = 300;

[D, x] = fourierDiff(N);
t = x*T/(2*pi);

Dmat = zeros(d*N);
Mmat = zeros(d*N);
for i = 1:N
    A = sysModel(t(i),zCoeff,a,b,M);
    for j = 1:d
       for k = 1:d
            Mmat(i+(j-1)*N,(k-1)*N+i) = A(j,k); 
       end
       Dmat((j-1)*N+1:j*N,(j-1)*N+1:j*N) = D*(2*pi/T);
    end
end

tic
[eigVec,eigE] = eig(Dmat-Mmat);
eigE = -1*diag(eigE);
eigM = exp(eigE*T);
t1 = toc;

% estimate FTM
tic
FTM = get_FTM(t,zCoeff,a,b,M,0);
FM = eig(FTM);
t2 = toc

FE = log(FM)/T;
display(FM)

save('eigData','eigM','FM','N','d','T','M')

rmpath('./getSoln_v2');

%%

figure
hold on
thet = linspace(0,2*pi,100);
plot(cos(thet),sin(thet),'-k','LineWidth',1);
axis equal
title_str = horzcat('Linear Wind , ','N = ',num2str(N));
title({'Floquet Multipliers',title_str});

tolVal = min(abs(FM))/5; 

idx = cell(d);
ix = [];
for i = 1:d
    idx{i} = find(abs(FM(i)-eigM)<tolVal);
    ix = [ix;idx{i}];
end
ix = unique(ix);
ixx = 1:length(eigM);
ixx(ix) = [];  
if isempty(ixx)
    p0 = plot(0,'m','LineWidth',1);
    p1 = plot(eigM(ix),'ob','LineWidth',1);
elseif isempty(ix)
    p0 = plot(eigM(ixx),'sm','LineWidth',1);
    p1 = plot(0,'b','LineWidth',1);
else
    p0 = plot(eigM(ixx),'sm','LineWidth',1);
    p1 = plot(eigM(ix),'ob','LineWidth',1);
end
p2 = plot(real(FM),imag(FM),'*r','LineWidth',1);
legend([p1,p0,p2],{horzcat('Accurate FM from spectral ','(',num2str(length(ix)),')'),horzcat('Spurious FM from spectral ','(',num2str(length(ixx)),')'),'FM from time-evolve'});

function A = sysModel(t,zCoeff,a,~,M)
    ff = fBasis(t,1,M);
    x = ff(1,:)*zCoeff(1,:)';
    y = ff(1,:)*zCoeff(2,:)';
    A = [-1+2*x*y, a+x^2;
           -2*x*y, -(a+x^2)];
end

% to calculate FTM
function FTM = get_FTM(t,zCoeff,a,b,M,flg)
global d
    switch flg
        case 0 % time-evolve method
            Phi=zeros(d);
            Phi0=eye(d);
            options = odeset('RelTol',1e-9,'AbsTol',1e-9);
            linDervs = @(t,X,zCoeff,a,b,M) sysModel(t,zCoeff,a,b,M)*X;
            for i=1:d
                [~,X]=ode15s(@(t,X) linDervs(t,X,zCoeff,a,b,M),[0,t(end)],Phi0(:,i),options);
                Phi(:,i)=X(end,:)';
            end
            FTM = Phi;
        case 1 % Friedmann's method
            h = t(2) - t(1);
            FTM = eye(d);
            for i = 1:length(t)-2
                K = friedmann_K(h, (t(end) - i*h),zCoeff,a,b,M);
                FTM = FTM*K;
            end
    end
 
end

function K = friedmann_K(h,t,zCoeff,a,b,M)
global d
    A_psi = sysModel(t,zCoeff,a,b,M);
    A_psi_h_by_2 = sysModel(t+0.5*h,zCoeff,a,b,M);
    A_psi_h = sysModel(t+h,zCoeff,a,b,M);
    E = A_psi_h_by_2*(eye(d) + 0.5*h*A_psi);
    F = A_psi_h_by_2*(eye(d) + (-0.5 + 2^(-0.5))*h*A_psi + (1 - 2^(-0.5))*h*E);
    G = A_psi_h*(eye(d) - h*(2^(-0.5))*E + (1 + 2^(-0.5))*h*F);       
    K = eye(d) + (h/6)*(A_psi + 2*(1 - 2^(-0.5))*E + 2*(1 + 2^(-0.5))*F + G);
end