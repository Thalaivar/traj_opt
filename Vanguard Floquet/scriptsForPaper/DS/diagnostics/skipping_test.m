addpath('..')
addpath('../lib')

clearvars
close all
load('..\full_state_discrete\solutions\EE50.mat')

type = 'full-state';
trajShape = 'eight';
%% make state and control matrices
if strcmp(type, 'full-state')
    x = zeros(N, 4);
    u = zeros(N, 2);
    for i = 1:8
        j = (i-1)*N + 1;
        if i <= 6
            x(:,i) = sol(j:j+N-1,1);
        else
            u(:,i-6) = sol(j:j+N-1,1);
        end
    end
    T = sol(8*N+1,1); VR = sol(8*N+2,1);
    [~,t] = fourierdiff(N);
    t = T*t/(2*pi);
    if strcmp(trajShape, 'circle')
        chiLinearTerm = sol(8*N+3);
        for i = 1:N
            x(i,2) = x(i,2) + chiLinearTerm*t(i);
        end
    end
end
traj_data.T = T; traj_data.VR = VR;
traj_data.state = x; traj_data.control = u; traj_data.p = p;
traj_data.t = t;
%% spectral floquet
global d;
n_spectral = N;
d = 6;
[D, x] = fourierdiff(n_spectral);
t = x*T/(2*pi);
Dmat = zeros(d*n_spectral);
Mmat = zeros(d*n_spectral);
for i = 1:n_spectral
    A = sysModel(t(i), traj_data);
    for j = 1:d
       for k = 1:d
            Mmat(i+(j-1)*n_spectral,(k-1)*n_spectral+i) = A(j,k); 
       end
       Dmat((j-1)*n_spectral+1:j*n_spectral,(j-1)*n_spectral+1:j*n_spectral) = D*(2*pi/T);
    end
end

[eigVec,eigE] = eig(Dmat-Mmat);
residVal = norm((Dmat-Mmat)*eigVec-eigVec*eigE);
% display(residVal)
eigE = -1*diag(eigE);
eigM = exp(eigE*T);
FTM = get_FTM(t, traj_data, 0);
FM = eig(FTM);
FE = log(FM)/T;

%% identification
eigE = eigE(find(abs(imag(eigE)) <= 0.25*N*2*pi/T));
sortedEigE = sort(eigE, 'ComparisonMethod', 'real');
scatter(real(sortedEigE), imag(sortedEigE), 'xm')

% identify vertical lines
verticalLines = []; j = 1; vLineTol = 1e-3; imCondTol = 1e-3;
eigGroups = {}; FE = [];
while(j+2 <= length(sortedEigE))
    j1 = sortedEigE(j); j2 = sortedEigE(j+1); j3 = sortedEigE(j+2);
    isVerticalLine = false;
    if abs(real(j1) - real(j2)) < vLineTol
        if abs(real(j2) - real(j3)) < vLineTol
            verticalLines(end+1) = real(j1);
            isVerticalLine = true;
        end
    end
    
    isFELine = false;
    if(isVerticalLine)
        imagCheck1 = mod(imag(j1) - imag(j2), 2*pi/T);
        imagCheck2 = mod(imag(j2) - imag(j3), 2*pi/T);
        imagCheck3 = mod(imag(j3) - imag(j1), 2*pi/T);
        if(imagCheck1 < imCondTol)
            isFELine = true;
        elseif(imagCheck2 < imCondTol)
            isFELine = true;
        elseif(imagCheck3 < imCondTol)
            isFELine = true;
        end
    end
    
    if(isFELine && isVerticalLine)
        currentEigGroup = [];
        verticalLineSize = 0;
        while(abs(real(sortedEigE(j)) - real(j1)) < vLineTol)
            currentEigGroup(end+1) = sortedEigE(j);
            verticalLineSize = verticalLineSize + 1; 
            j = j + 1;
            if j > length(sortedEigE)
                break
            end
        end
        eigGroups{end+1} = currentEigGroup;
        if verticalLineSize < 0.75*N
            FE(end+1) = real(j1);
        else
            if verticalLineSize > 1.25*N
                FE(end+1) = real(j1);
                FE(end+1) = real(j1);
                FE(end+1) = real(j1);
            else
                FE(end+1) = real(j1);
                FE(end+1) = real(j1);
            end
        end
    else
        j = j + 1;
    end
end

rmpath('..')
rmpath('../lib')

% Conclusions:
% 1. Two eigenvalues with real parts same uptil the 4th decimal place are
% considered to be belonging to the same vertical line.
% 2. 
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
    Z = [V, chi, gamma, x, y, z];
    U = [CL, mu, 0, traj_data.VR];
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
        case 1 % Friedman_spectral's method
            h = t(2) - t(1);
            FTM = eye(d);
            for i = 1:length(t)-2
                K = friedman_spectral_K(h, (t(end) - i*h));
                FTM = FTM*K;
            end
    end
 
end

function K = friedman_spectral_K(h,t)
global d
    A_psi = sysModel(t);
    A_psi_h_by_2 = sysModel(t+0.5*h);
    A_psi_h = sysModel(t+h);
    E = A_psi_h_by_2*(eye(d) + 0.5*h*A_psi);
    F = A_psi_h_by_2*(eye(d) + (-0.5 + 2^(-0.5))*h*A_psi + (1 - 2^(-0.5))*h*E);
    G = A_psi_h*(eye(d) - h*(2^(-0.5))*E + (1 + 2^(-0.5))*h*F);       
    K = eye(d) + (h/6)*(A_psi + 2*(1 - 2^(-0.5))*E + 2*(1 + 2^(-0.5))*F + G);
end