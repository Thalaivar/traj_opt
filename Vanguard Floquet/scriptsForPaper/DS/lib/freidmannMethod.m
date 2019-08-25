% trajData is struct with fields:
%   - T (time period)
%   - fourierGrid (the fourier grid (0, 2*pi] for the trajectory
%   - X (Nx6 state matrix)
%   - U (Nx3 control matrix)
%   - type ('diff-flat' for flatness or 'full-state')
%   - VR (reference wind speed)
function [FTM,FE] = freidmannMethod(trajData, freidmannGrid)
    t = trajData.T*freidmannGrid/(2*pi);
    d = 4;
    h = t(2) - t(1);
    FTM = eye(d);
    for i = 1:length(t)-2
        K = friedman_spectral_K(h, (t(end) - i*h), d, trajData);
        FTM = FTM*K;
    end
    FM = eig(FTM);
    FE = log(FM)/trajData.T;
end

function K = friedman_spectral_K(h,t,d,trajData)
    A_psi = sysModel(t, trajData);
    A_psi_h_by_2 = sysModel(t+0.5*h,trajData);
    A_psi_h = sysModel(t+h, trajData);
    E = A_psi_h_by_2*(eye(d) + 0.5*h*A_psi);
    F = A_psi_h_by_2*(eye(d) + (-0.5 + 2^(-0.5))*h*A_psi + (1 - 2^(-0.5))*h*E);
    G = A_psi_h*(eye(d) - h*(2^(-0.5))*E + (1 + 2^(-0.5))*h*F);       
    K = eye(d) + (h/6)*(A_psi + 2*(1 - 2^(-0.5))*E + 2*(1 + 2^(-0.5))*F + G);
end

function A = sysModel(t, trajData)
    prm.m = 4.5;
    prm.S = 0.473;
    prm.CD0 = 0.0173;
    prm.CD1 = -0.0337;
    prm.CD2 = 0.0517;
    prm.p_exp = trajData.p;
    % get state and control vector for jacobian
    tt  = trajData.T*trajData.fourierGrid/(2*pi);
    V     = interp_sinc(tt, trajData.X(:,1), t);
    chi   = interp_sinc(tt, trajData.X(:,2), t);
    gamma = interp_sinc(tt, trajData.X(:,3), t);
    x     = interp_sinc(tt, trajData.X(:,4), t);
    y     = interp_sinc(tt, trajData.X(:,5), t);
    z     = interp_sinc(tt, trajData.X(:,6), t);
    CL = interp_sinc(tt, trajData.U(:,1), t);
    mu = interp_sinc(tt, trajData.U(:,2), t);
    if strcmp(trajData.type, 'diff-flat')
        CT = interp_sinc(tt, trajData.U(:,3), t);
    else
        CT = 0;
    end
    Z = [V, chi, gamma, x, y, z];
    U = [CL, mu, CT, trajData.VR];
    % evaluate jacobian
    A = JacEval(Z, U, prm);
    A = [A(1:3,1:3), A(1:3,6); A(6,1:3), A(6,6)];
end
