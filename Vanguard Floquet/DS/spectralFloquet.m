addpath('lib\')

clearvars
close all
load('EET_200.mat')
x = zeros(N, 4);
u = zeros(N, 3);
for i = 1:9
    j = (i-1)*N + 1;
    if i <= 6
        x(:,i) = sol(j:j+N-1,1);
    else
        u(:,i-6) = sol(j:j+N-1,1);
    end
end
T = sol(9*N+1,1); VR = sol(9*N+2,1);
traj_data.T = T; traj_data.VR = VR;
traj_data.state = x; traj_data.control = u; traj_data.p = p;
%% spectral floquet
global d;
NN = 150;
d = 6;
[D, x] = fourierdiff(NN);
t = x*T/(2*pi);
traj_data.t = t;
Dmat = zeros(d*NN);
Mmat = zeros(d*NN);
for i = 1:NN
    A = sysModel(t(i), traj_data);
    for j = 1:d
       for k = 1:d
            Mmat(i+(j-1)*NN,(k-1)*NN+i) = A(j,k); 
       end
       Dmat((j-1)*NN+1:j*NN,(j-1)*NN+1:j*NN) = D*(2*pi/T);
    end
end

[eigVec,eigE] = eig(Dmat-Mmat);
residVal = norm((Dmat-Mmat)*eigVec-eigVec*eigE);
display(residVal)
eigE = -1*diag(eigE);
eigM = exp(eigE*T);

FTM = get_FTM(t, traj_data, 0);
FM = eig(FTM);
FE = log(FM)/T;

figure
hold on
p1 = plot(real(eigE),imag(eigE),'xm');
p2 = plot(real(FE),imag(FE),'or','LineWidth',1,'MarkerSize',10);
limSet = 1.2*max(abs(real(eigE)));
xlim([-limSet,limSet]);
if p == 1
    traj_type = 'linear';
else
    traj_type = 'exponential';
end
title_str1 = ['Dynamic soaring in ', traj_type, ' profile. '];
title_str2 = ['N = ' num2str(NN)];
title([title_str1 title_str2]);
grid minor
legend([p1,p2],{'Spectral FE','Time-evolved FE'});

estimated_FE = identify_floquet(eigE, N, traj_data.T);
scatter(real(estimated_FE), imag(estimated_FE), 'sb', 'LineWidth', 1.25, 'DisplayName', 'Estimated FE')
%% supporting functions
    function A = sysModel(t, traj_data)
    prm.m = 4.5;
    prm.S = 0.473;
    prm.CD0 = 0.0173;
    prm.CD1 = -0.0337;
    prm.CD2 = 0.0517;
    prm.p_exp = traj_data.p;
    % get state and control vector for jacobian
    V     = interp_sinc(traj_data.t, traj_data.state(:,1), t);
    chi   = interp_sinc(traj_data.t, traj_data.state(:,2), t);
    gamma = interp_sinc(traj_data.t, traj_data.state(:,3), t);
    x     = interp_sinc(traj_data.t, traj_data.state(:,4), t);
    y     = interp_sinc(traj_data.t, traj_data.state(:,5), t);
    z     = interp_sinc(traj_data.t, traj_data.state(:,6), t);
    CL = interp_sinc(traj_data.t, traj_data.control(:,1), t);
    mu = interp_sinc(traj_data.t, traj_data.control(:,2), t);
    CT = interp_sinc(traj_data.t, traj_data.control(:,3), t);
    Z = [V, chi, gamma, x, y, z];
    U = [CL, mu, CT, traj_data.VR];
    % evaluate jacobian
    A = JacEval(Z, U, prm);
end

function FTM = get_FTM(t,traj_data, flg)
global d
    switch flg
        case 0 % time-evolve method
            Phi=zeros(d);
            Phi0=eye(d);
            options = odeset('RelTol',1e-9,'AbsTol',1e-9);
            linDervs = @(t,X) sysModel(t, traj_data)*X;
            for i=1:d
                [~,X] = ode15s(@(t,X) linDervs(t,X),[0,t(end)],Phi0(:,i),options);
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

function K = friedmann_K(h,t)
global d
    A_psi = sysModel(t);
    A_psi_h_by_2 = sysModel(t+0.5*h);
    A_psi_h = sysModel(t+h);
    E = A_psi_h_by_2*(eye(d) + 0.5*h*A_psi);
    F = A_psi_h_by_2*(eye(d) + (-0.5 + 2^(-0.5))*h*A_psi + (1 - 2^(-0.5))*h*E);
    G = A_psi_h*(eye(d) - h*(2^(-0.5))*E + (1 + 2^(-0.5))*h*F);       
    K = eye(d) + (h/6)*(A_psi + 2*(1 - 2^(-0.5))*E + 2*(1 + 2^(-0.5))*F + G);
end