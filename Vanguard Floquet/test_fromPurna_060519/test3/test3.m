% 09/04/19
clear all

addpath('./JacEval');
addpath('./AllData');
load AllData

global splStates splControls
VR = data_traj.VR;
splStates = spline(data_traj.t',data_traj.y');
splControls = spline(data_traj.t',data_traj.u');

global d
d = 6;
N = 50;

[D, x] = fourierDiff(N);
T = data_traj.tf;
t = x*T/(2*pi);

Dmat = zeros(d*N);
M = zeros(d*N);
for i = 1:N
    A = sysModel(t(i),prm,VR);
    for j = 1:d
       for k = 1:d
            M(i+(j-1)*N,(k-1)*N+i) = A(j,k); 
       end
       Dmat((j-1)*N+1:j*N,(j-1)*N+1:j*N) = D*(2*pi/T);
    end
end

% Dmat = sparse(Dmat);
% M = sparse(M);

[eigVec,eigE] = eig(Dmat-M);
eigE = -1*diag(eigE);
eigM = exp(eigE*T);

% estimate FTM
FTM = get_FTM(t,prm,VR,0);
FM = eig(FTM);
FE = log(FM)/T;
 
rmpath('./JacEval');
rmpath('./AllData')

% save('eigData','eigM','FM','N','d','T','M')

%%

figure
hold on
thet = linspace(0,2*pi,100);
plot(cos(thet),sin(thet),'-k','LineWidth',1);
axis equal
title_str = horzcat('Linear Wind , ','N = ',num2str(N));
title({'Floquet Multipliers',title_str});

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
    
function A = sysModel(t,prm,VR)
global splStates splControls
X = ppval(splStates,t);
U = ppval(splControls,t);
U = [U;VR];
A = JacEval(X,U,prm);
end

% to calculate FTM
function FTM = get_FTM(t,prm,VR,flg)
global d
    switch flg
        case 0 % time-evolve method
            Phi=zeros(d);
            Phi0=eye(d);
            options = odeset('RelTol',1e-9,'AbsTol',1e-9);
            linDervs = @(t,X,prm,VR) sysModel(t,prm,VR)*X;
            for i=1:d
                [~,X]=ode15s(@(t,X) linDervs(t,X,prm,VR),[0,t(end)],Phi0(:,i),options);
                Phi(:,i)=X(end,:)';
            end
            FTM = Phi;
        case 1 % Friedmann's method
            h = t(2) - t(1);
            FTM = eye(d);
            for i = 1:length(t)-2
                K = friedmann_K(h, (t(end) - i*h),prm,VR);
                FTM = FTM*K;
            end
    end
 
end

function K = friedmann_K(h,t,prm,VR)
global d
    A_psi = sysModel(t,prm,VR);
    A_psi_h_by_2 = sysModel(t+0.5*h,prm,VR);
    A_psi_h = sysModel(t+h,prm,VR);
    E = A_psi_h_by_2*(eye(d) + 0.5*h*A_psi);
    F = A_psi_h_by_2*(eye(d) + (-0.5 + 2^(-0.5))*h*A_psi + (1 - 2^(-0.5))*h*E);
    G = A_psi_h*(eye(d) - h*(2^(-0.5))*E + (1 + 2^(-0.5))*h*F);       
    K = eye(d) + (h/6)*(A_psi + 2*(1 - 2^(-0.5))*E + 2*(1 + 2^(-0.5))*F + G);
end